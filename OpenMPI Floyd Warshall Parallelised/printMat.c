/**
 * Prints Binary matrix
 */
#include <stdio.h>
#include <stdlib.h>

int fileRead(char* file) {
    FILE* fd = fopen(file, "rb");
    int numV;

    if (fd == NULL) return -1;

    fread(&numV, sizeof(int), 1, fd);
    int matSize = (numV) * (numV);
    printf("%d\n", numV);

    int temp2 = 0;
    for (int i = 0; i < matSize; i++) {
        fread(&temp2, sizeof(int), 1, fd);
        printf("%d\n", temp2);
    }

    fclose(fd);
    return 0;
}

int main(int argc, char** argv) {
    fileRead(argv[1]);
    return 0;
}
