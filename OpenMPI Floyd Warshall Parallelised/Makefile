# Makefile to make project2.c

CC=mpicc
CFLAGS=-Wall -Werror -pedantic -std=c99
EXE=./bin/project2.exe

project2: project2.c # compile project
	$(CC) $(CFLAGS) -o $(EXE) project2.c

run: $(EXE) # Run executable with mpirun
	mpirun -n 4 $(EXE) $< $(ARGS)
