#! /usr/bin/bash

cd tlb-nh-optimized
sh test_numa.sh
cd no_prefetch
sh test_nvm.sh
sh test_o2.sh
#sh test.sh
#cd ../../
#sh test_nvm.sh
#sh test_pmm.sh
