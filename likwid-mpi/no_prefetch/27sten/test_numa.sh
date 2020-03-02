#!/usr/bin/bash

make
export OMP_NUM_THREADS=1
likwid-mpirun -np 48 -g L3CACHE -m -C 0-47	./numatest_mpi 4294967296 &> numa_ss_np4_48_l3
likwid-mpirun -np 48 -g CACHES -m -C 0-47	./numatest_mpi 4294967296 &> numa_ss_np4_48_cb
likwid-mpirun -np 48 -g TLB_DATA -m	 -C 0-47 ./numatest_mpi 4294967296 &> numa_ss_np4_48_td
