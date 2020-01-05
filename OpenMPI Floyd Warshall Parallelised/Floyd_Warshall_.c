/**
 * CITS3402 Project 2: Shortest Paths
 * Parallelised Floyd Warshall
 *
 * Written by Nicholas Cannon (22241579), Paul O'Sullivan (21492318)
 */
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define LOG_PATH "logs/"
#define COLLECT_MSG 42
#define INFINITY 100000

/**
 * Reads input file, adding the contents to a 1 dimensional array.
 * Returns non-zero on failure.
 */
int fileRead(char* file, int** distance, int* numV) {
    FILE* fd = fopen(file, "rb");

    if (fd == NULL) return -1;

    fread(numV, sizeof(int), 1, fd);
    int matSize = (*numV) * (*numV);

    *distance = malloc(matSize * sizeof(int));

    int temp2 = 0;
    for (int i = 0; i < matSize; i++) {
        fread(&temp2, sizeof(int), 1, fd);
        (*distance)[i] = temp2;
    }

    fclose(fd);
    return 0;
}

/**
 * Writes distance matrix to file. Format is LOG_PATH + <DDMMYYYY>_<HHMM>.out
 * Returns non-zero on failure.
 */
int writeFile(int* distance, int numV) {
    // Format time string
    time_t timer;
    time(&timer);
    struct tm* tm_info;
    tm_info = localtime(&timer);
    char timestr[14];
    strftime(timestr, 14, "%d%m%Y_%H%M", tm_info);

    // form log path
    char* logPath =
        malloc(sizeof(char) * (strlen(LOG_PATH) + strlen(timestr) + 5));
    strcpy(logPath, LOG_PATH);
    strcat(logPath, timestr);
    strcat(logPath, ".out");

    FILE* fd = fopen(logPath, "wb");
    if (fd == NULL) return -1;

    fwrite(&numV, sizeof(int), 1, fd);  // write size of matrix
    for (int i = 0; i < numV * numV; i++) {
        fwrite(&distance[i], sizeof(int), 1, fd);
    }

    fclose(fd);
    free(logPath);
    return 0;
}

/**
 * Returns the number of rows owned by process.
 *  rank - rank of the process
 *  numV - number of vertices in graph
 *  numP - number of processes
 */
int getNumRows(int rank, int numV, int numP) {
    return (((rank + 1) * numV) / numP) - ((rank * numV) / numP);
}

/**
 * Returns the rank of the process that owns given row.
 *  row  -  row in matrix
 *  numV -  number of vertices in graph
 *  size -  number of processes in communicator
 */
int getOwner(int row, int numV, int size) {
    return (size * (row + 1) - 1) / numV;
}

/**
 * Parallel implementation of Floyd Warshall all pairs shortest path.
 *  rank        -   rank of the process running
 *  size        -   total number of processes
 *  numV        -   number of vertices in the graph
 *  localRows   -   processes local partition of the matrix
 */
void parallelFW(int rank, int size, int numV, int* localRows) {
    int* temp = malloc(sizeof(int) * numV);
    int owner;

    for (int row = 0; row < numV; row++) {
        owner = getOwner(row, numV, size);
        // check if this process owns this row
        if (rank == owner) {
            // broadcast this row to all other processes
            int subIndex = row - ((rank * numV) / size);  // relative index
            for (int i = 0; i < numV; i++) {
                temp[i] = localRows[(numV * subIndex) + i];
            }
        }
        // have owner broadcast row to other processes
        MPI_Bcast(temp, numV, MPI_INT, owner, MPI_COMM_WORLD);
        int rowCount = getNumRows(rank, numV, size);
        for (int lRow = 0; lRow < rowCount; lRow++) {
            for (int element = 0; element < numV; element++) {
                int newLength = localRows[(lRow * numV) + row] + temp[element];
                if (newLength < localRows[(lRow * numV) + element]) {
                    localRows[(lRow * numV) + element] = newLength;
                }
            }
        }
    }
    free(temp);
}

int main(int argc, char** argv) {
    int rank, size, status, numV, rows, msg = COLLECT_MSG;
    int* distance = NULL;
    int* localRows = NULL;
    double start, end;

    // Quick argument error check
    if (argc < 3) {
        printf("Please supply a graph with -f option\n");
        return (EXIT_FAILURE);
    }

    // Initialize our MPI environment and get process specific data
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Parse command line arguments
        char* filename = NULL;
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "-f") == 0) {
                if (argv[i + 1]) {
                    filename = argv[i + 1];
                }
            }
        }

        // read in the whole graph
        status = fileRead(filename, &distance, &numV);
        if (status != 0) {
            // error occurred reading graph
            printf("Error reading graph: %s\n", filename);
            numV = -1;
        }

        // remove cyclic paths and add "infinity" for all non-existent paths
        for (int i = 0; i < numV; i++) {
            for (int j = 0; j < numV; j++) {
                if (distance[(i * numV) + j] == 0) {
                    // non-existent path
                    distance[(i * numV) + j] = INFINITY;
                }
                if (i == j) {
                    distance[(i * numV) + j] = 0;
                }
            }
        }

        // start wall time
        start = MPI_Wtime();
    }

    // send matrix size to every process
    MPI_Bcast(&numV, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (numV == -1) {
        // error occurred when reading in matrix -> abort all processes
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // allocated memory for local partition of matrix
    rows = getNumRows(rank, numV, size);
    localRows = malloc(sizeof(int) * rows * numV);

    // Partition data
    if (rank == 0) {
        // give master data before other processes
        int sent = rows * numV;
        memcpy(localRows, distance, sizeof(int) * sent);

        // Partition matrix and send to other processes
        for (int r = 1; r < size; r++) {
            int elementsToSend = getNumRows(r, numV, size) * numV;
            MPI_Send((distance + sent), elementsToSend, MPI_INT, r, 0,
                     MPI_COMM_WORLD);
            sent += elementsToSend;
        }
    } else {
        // recieve partition of matrix from master
        MPI_Recv(localRows, rows * numV, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }

    // Process local rows here
    parallelFW(rank, size, numV, localRows);

    // Gather and log shortest paths to file
    if (rank == 0) {
        // load master process data into distance array
        for (int i = 0; i < rows * numV; i++) {
            distance[i] = localRows[i];
        }

        // gather rows from other processes
        int recieved = rows * numV;
        for (int r = 1; r < size; r++) {
            // let process r know we want their data and recieve it
            MPI_Send(&msg, 1, MPI_INT, r, 0, MPI_COMM_WORLD);
            MPI_Recv((distance + recieved), getNumRows(r, numV, size) * numV,
                     MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recieved += getNumRows(r, numV, size) * numV;
        }

        end = MPI_Wtime();
        printf("time taken = %f sec\n", end - start);

        // log the matrix
        status = writeFile(distance, numV);
        if (status != 0) printf("Error logging output\n");
        printf("Written log file\n");
    } else {
        // wait until master asks for our data then send it
        int ans;
        MPI_Recv(&ans, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(localRows, rows * numV, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    // clean up, terminate the MPI environment and return
    free(distance);
    free(localRows);
    MPI_Finalize();
    return 0;
}
