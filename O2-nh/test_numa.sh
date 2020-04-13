#!/usr/bin/bash

make
export OMP_NUM_THREADS=1
mpirun --bind-to-core -np 24	./numatest_mpi 4294967296 &> numa_ss_pre4_24
mpirun --bind-to-core -np 48	./numatest_mpi 4294967296 &> numa_ss_pre4_48
