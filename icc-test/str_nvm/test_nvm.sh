#!/usr/bin/bash

make
export OMP_NUM_THREADS=1
likwid-mpirun -np 48 -g PMM -m 	./numatest_mpi 4294967296 &> numa_ss_np4_48_pmm
