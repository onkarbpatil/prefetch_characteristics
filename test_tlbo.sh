#! /usr/bin/bash

cd tlb-optimized
sh test_nvm.sh
cd no_prefetch
sh test_numa.sh
sh test_nvm.sh
sh test_o2.sh
#sh test.sh
#cd ../../
#sh test_nvm.sh
#sh test_pmm.sh
