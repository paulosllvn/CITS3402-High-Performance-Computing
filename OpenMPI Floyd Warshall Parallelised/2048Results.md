# Speed up (wall time) of 2048 graph

## 1 process per node

```
mpirun --hostfile hostfiles/host1 -N 1 project2 -f graphs/2048.in
time taken = 55.189667 sec

mpirun --hostfile hostfiles/host2 -N 1 project2 -f graphs/2048.in
time taken = 27.696308 sec

mpirun --hostfile hostfiles/host4 -N 1 project2 -f graphs/2048.in
time taken = 14.358642 sec

mpirun --hostfile hostfiles/host8 -N 1 project2 -f graphs/2048.in
time taken = 7.699995 sec

mpirun --hostfile hostfiles/host16 -N 1 project2 -f graphs/2048.in
time taken = 4.560019 sec

```
