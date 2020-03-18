#! /usr/bin/bash

cd tlb-prefetch
sh test_numa.sh
cd no_prefetch
sh test_nvm.sh
#sh test.sh
#cd ../../
#sh test_nvm.sh
#sh test_pmm.sh
