#!/usr/bin/bash

rm numatest_mpi
make
export OMP_NUM_THREADS=1
likwid-mpirun -np 48 -g L3CACHE -m 	./numatest_mpi 4294967296 &> numa_ss_pre4_48_l3
likwid-mpirun -np 48 -g CACHES -m 	./numatest_mpi 4294967296 &> numa_ss_pre4_48_cb
likwid-mpirun -np 48 -g TLB_DATA -m  ./numatest_mpi 4294967296 &> numa_ss_pre4_48_td
