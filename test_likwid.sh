#! /usr/bin/bash

cd likwid-mpi
sh test_numa.sh
cd no_prefetch
sh test_numa.sh
sh test.sh
cd ../../
sh test_nvm.sh
