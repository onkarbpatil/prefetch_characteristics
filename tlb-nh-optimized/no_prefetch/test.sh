#! /usr/bin/bash

cd wr-only
sh test_numa.sh
cd ../1w4r
sh test_numa.sh
cd ../str
sh test_numa.sh
cd ../rand
sh test_numa.sh
cd ../l3cache
sh test_numa.sh
cd ../3sten
sh test_numa.sh
cd ../5sten
sh test_numa.sh
cd ../7sten
sh test_numa.sh
cd ../9sten
sh test_numa.sh
cd ../27sten
sh test_numa.sh