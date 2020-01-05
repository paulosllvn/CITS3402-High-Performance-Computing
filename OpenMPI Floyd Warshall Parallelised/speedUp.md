# Graph Speed up

## 4 nodes -> 1 process per node

```
mpirun --hostfile hostfiles/host4 -N 1 project2 -f graphs/512.in
time taken = 0.294059 sec

mpirun --hostfile hostfiles/host4 -N 1 project2 -f graphs/1024.in
time taken = 2.040283 sec

mpirun --hostfile hostfiles/host4 -N 1 project2 -f graphs/2048.in
time taken = 14.436905 sec

mpirun --hostfile hostfiles/host4 -N 1 project2 -f graphs/4096.in
time taken = 111.311171 sec
```

## 2 nodes -> 2 processes per node (4 processes total)

```
mpirun --hostfile hostfiles/host2 -N 2 project2 -f graphs/512.in
time taken = 0.292934 sec

mpirun --hostfile hostfiles/host2 -N 2 project2 -f graphs/1024.in
time taken = 2.088377 sec

mpirun --hostfile hostfiles/host2 -N 2 project2 -f graphs/2048.in
time taken = 14.428202 sec

mpirun --hostfile hostfiles/host2 -N 2 project2 -f graphs/4096.in
time taken = 115.155423 sec
```
