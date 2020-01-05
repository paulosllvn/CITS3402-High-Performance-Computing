/**
 * Sequential implementation of Floyd Warshall.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define LOG_PATH "logs/seq/"

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
    return 0;
}

/**
 * Performs a sequential floyd Warshall shortest path algorithm.
 */
void seqFloydWarshall(int* distance, int* sp, int numV) {
    for (int k = 0; k < numV; k++) {
        for (int i = 0; i < numV; i++) {
            for (int j = 0; j < numV; j++) {
                if (distance[(i * numV) + j] >
                    distance[(i * numV) + k] + distance[(k * numV) + j]) {
                    sp[(i * numV) + j] =
                        distance[(i * numV) + k] + distance[(k * numV) + j];
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    int* dist = NULL;
    char* filename = NULL;
    int status, numV;

    // get command line options
    if (argc < 3) {
        printf("Invalid usage\n");
        return EXIT_FAILURE;
    }

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0) {
            if (argv[i + 1]) {
                filename = argv[i + 1];
            }
        }
    }

    // read in matrix
    status = fileRead(filename, &dist, &numV);
    if (status != 0) {
        printf("Error reading matrix in %s\n", filename);
        return EXIT_FAILURE;
    }

    // remove cyclic paths
    for (int i = 0; i < numV; i++) {
        for (int j = 0; j < numV; j++) {
            if (dist[(i * numV) + j] == 0) {
                // non-existent path
                dist[(i * numV) + j] = 100000;
            }
            if (i == j) {
                dist[(i * numV) + j] = 0;
            }
        }
    }

    // calc shortest paths with FW
    int* sp = malloc(sizeof(int) * numV * numV);
    memcpy(sp, dist, numV * numV * sizeof(int));
    seqFloydWarshall(dist, sp, numV);

    for (int i = 0; i < numV; i++) {
        for (int j = 0; j < numV; j++) {
            printf("%d", sp[(i * numV) + j]);
        }
    }

    // write log file
    // status = writeFile(dist, numV);
    // if (status != 0) {
    //     printf("Error writing log file\n");
    //     return EXIT_FAILURE;
    // }

    free(dist);
    return 0;
}
