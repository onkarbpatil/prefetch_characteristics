!/usr/bin/bash

make
mpirun --bind-to-core -np 24	./numatest_mpi 4294967296 &> numa_ss_pre4_mm-sonata_24
mpirun --bind-to-core -np 48	./numatest_mpi 4294967296 &> numa_ss_pre4_mm-sonata_48
