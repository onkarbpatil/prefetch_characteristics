#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <sched.h>
#include <sys/time.h>
#include "numa.h"
#include <mpi.h>

#define BILLION  1000000000L;

int mem_types;
int max_node;
int numt;
int total_numa_nodes;
int * numa_node_ids;
struct bitmask * numa_nodes;
char ** mem_tech;
long double * means;
int * cluster_sizes;

void numatest(int argc, char ** argv, int rank, int procs, unsigned long bytes);
