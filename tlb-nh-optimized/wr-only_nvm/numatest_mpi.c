#include "numatest_mpi.h"

int mem_types;
int max_node;
int numt;
int total_numa_nodes = 0;
int * numa_node_ids;
struct bitmask * numa_nodes;



void numatest(int argc, char ** argv, int rank, int procs, unsigned long bytes){
	max_node = numa_max_node() + 1;
	int cpu_count = numa_num_possible_cpus();
	numa_node_ids = (int*)malloc(sizeof(int)*max_node);
	struct bitmask * numa_nodes = numa_get_membind();
	int i = 0;
        while(i < numa_nodes->size){
                if(numa_bitmask_isbitset(numa_nodes, i)){
                        numa_node_ids[total_numa_nodes] = i;
                        total_numa_nodes++;
                }
                i++;
        }
		total_numa_nodes++;
	unsigned long size = bytes/procs;
	int mbs = size/sizeof(double);
	int ldim = (int)cbrt((double)size/sizeof(double))-1;
	int n = ldim;
	int rs = 0;
	int z = 0;
	int bdist = 0;
	int sdist = 0;
	int rd_dist, wr_dist;
	int unroll = 0;
	int * rand_tab;
	int inner = 0;
	int ext_flg = 0;
	int pg_sz = 4096;
	rand_tab = (int*)malloc(mbs*sizeof(int));
	double *a, *b, *c, *d, *e, *f, *g, *h;
	double **aa, **bb, **cc;
	clock_t start, end;
	struct timespec begin, stop, obegin, ostop;
	srand(10725);
	if(argc == 0){
		printf("Enter memory technologies available in ascending order of speed. eg: GPU NVRAM DRAM HBM\n");
		return;
	}
	else{
		mem_types = argc;
		mem_tech = (char**)malloc(argc*sizeof(char*));
		int a;
		for(a = 0; a < argc; a++){
			mem_tech[a] = argv[a];
		}

	}
	for(i=0; i< size/sizeof(double); i++)
		rand_tab[i]=rand()%(size/sizeof(double));

  	i = 2;
//	while(i < total_numa_nodes){
	{
wr_dist = 4096/sizeof(double);
//while(wr_dist < 32768/sizeof(double)){
{
	rd_dist = 32/sizeof(double);
//		while(rd_dist < 32768/sizeof(double)){
	{
				unroll = 64;
//				while(unroll < 128){
				{
	// Dynamically allocate the three arrays using "posix_memalign()"
		int iters = 0;
		int stride;
		long double wr_only_avg=0.0;
		long double owfr_avg=0.0;
		long double str_avg=0.0;
		long double rand_avg = 0.0;
		long double l2cache_avg = 0.0;
		long double t_sten_avg = 0.0;
		long double f_sten_avg = 0.0;
		long double s_sten_avg = 0.0;
		long double n_sten_avg = 0.0;
		long double t7_sten_avg = 0.0;
		long double wr_only_t = 0.0;
		long double owfr_t = 0.0;
		long double l2cache_t = 0.0;
		long double str_t = 0.0;
		long double rand_t = 0.0;
		long double t_sten_t = 0.0;
		long double f_sten_t = 0.0;
		long double s_sten_t = 0.0;
		long double n_sten_t = 0.0;
		long double t7_sten_t = 0.0;
		long double accum;

//		for( iters = 0; iters < 10; iters++)
		{
			int j = 0;
			int k = 0;
			int l = 0;
/*			if(i == (total_numa_nodes-1)){
				a = (double*)numa_alloc_onnode(size, numa_node_ids[0]);
				b = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				c = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				d = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
				e = (double*)numa_alloc_onnode(size, numa_node_ids[2]);
			}else{*/
				a = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				b = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				c = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				d = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
				e = (double*)numa_alloc_onnode(size, numa_node_ids[i]);
		//	}
			long double empty=0.0;
			long double empty2=0.0;

			MPI_Barrier(MPI_COMM_WORLD);
			LIKWID_MARKER_START("reg");
			clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
			for(j = 0;j < wr_dist;j++){
					__builtin_prefetch (&a[j], 1, 0);
					__builtin_prefetch (&b[j], 1, 0);
					__builtin_prefetch (&c[j], 1, 0);
					__builtin_prefetch (&d[j], 1, 0);
					__builtin_prefetch (&e[j], 1, 0);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(j = 0;j < ((size/sizeof(double)) - wr_dist);j+=unroll){
					if((wr_dist < pg_sz)||(j%(pg_sz/sizeof(double)) == 0)){
					__builtin_prefetch (&a[j+wr_dist], 1, 0);
					__builtin_prefetch (&b[j+wr_dist], 1, 0);
					__builtin_prefetch (&c[j+wr_dist], 1, 0);
					__builtin_prefetch (&d[j+wr_dist], 1, 0);
					__builtin_prefetch (&e[j+wr_dist], 1, 0);
					}
					for(inner = 0;inner<unroll;inner++){
						a[j+inner] = 1.0;
						b[j+inner] = 2.0;
						c[j+inner] = 3.0;
						d[j+inner] = 4.0;
						e[j+inner] = 5.0;
					}
			}
			for(j = ((size/sizeof(double)) - wr_dist); j < (size/sizeof(double)); j++){
				a[j] = 1.0;
				b[j] = 2.0;
				c[j] = 3.0;
				d[j] = 4.0;
				e[j] = 5.0;
			}
			LIKWID_MARKER_STOP("reg");
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			wr_only_t += accum;
			wr_only_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
/*
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
			if(wr_dist >= rd_dist){
				for(j =0; j < wr_dist - rd_dist; j++){
					__builtin_prefetch (&a[j], 1, 0);
				}
				for(j =(wr_dist - rd_dist); j < wr_dist; j++){
					__builtin_prefetch (&a[j], 1, 0);
					__builtin_prefetch (&b[j-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&c[j-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&d[j-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&e[j-(wr_dist - rd_dist)], 0, 0);
				}
			}else{
				for(j =0; j < rd_dist - wr_dist; j++){
					__builtin_prefetch (&b[j], 0, 0);
					__builtin_prefetch (&c[j], 0, 0);
					__builtin_prefetch (&d[j], 0, 0);
					__builtin_prefetch (&e[j], 0, 0);
				}
				for(j =(rd_dist - wr_dist); j < rd_dist; j++){
					__builtin_prefetch (&a[j-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[j], 0, 0);
					__builtin_prefetch (&c[j], 0, 0);
					__builtin_prefetch (&d[j], 0, 0);
					__builtin_prefetch (&e[j], 0, 0);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(j =0; j < ((size/sizeof(double))-bdist); j+=unroll){
					if((wr_dist < pg_sz)||(j%(pg_sz/sizeof(double)) == 0)){
					__builtin_prefetch (&a[j+wr_dist], 1, 0);
					}
					
					if((rd_dist < pg_sz)||(j%(pg_sz/sizeof(double)) == 0)){
					__builtin_prefetch (&b[j+rd_dist], 0, 0);
					__builtin_prefetch (&c[j+rd_dist], 0, 0);
					__builtin_prefetch (&d[j+rd_dist], 0, 0);
					__builtin_prefetch (&e[j+rd_dist], 0, 0);
					}
						for(inner=0;inner<unroll;inner++)
                            a[j+inner] = c[j+inner] + d[j+inner] +e[j+inner] + b[j+inner];
            }
			if(wr_dist >= rd_dist){
				for(j = ((size/sizeof(double))-wr_dist);j <((size/sizeof(double))-rd_dist);j+=unroll){
					if((rd_dist < pg_sz)||(j%(pg_sz/sizeof(double)) == 0)){
						__builtin_prefetch (&b[j+rd_dist], 0, 0);
						__builtin_prefetch (&c[j+rd_dist], 0, 0);
						__builtin_prefetch (&d[j+rd_dist], 0, 0);
						__builtin_prefetch (&e[j+rd_dist], 0, 0);
					}
						for(inner=0;inner<unroll;inner++)
						a[j+inner] = c[j+inner] + d[j+inner] + e[j+inner] + b[j+inner];
				}
				for(j = ((size/sizeof(double))-rd_dist);j <(size/sizeof(double));j++){
						a[j] = c[j] + d[j] + e[j] + b[j];
				}
			}else{
				for(j = ((size/sizeof(double))-rd_dist);j <((size/sizeof(double))-wr_dist);j+=unroll){
					if((wr_dist < pg_sz)||(j%(pg_sz/sizeof(double)) == 0)){
						__builtin_prefetch (&a[j+wr_dist], 1, 0);
					}
						for(inner=0;inner<unroll;inner++)
						a[j+inner] = c[j+inner] + d[j+inner] + e[j+inner] + b[j+inner];
				}
				for(j = ((size/sizeof(double))-wr_dist);j <(size/sizeof(double));j++){
						a[j] = c[j] + d[j] + e[j] + b[j];
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			owfr_t += accum;
			owfr_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}

			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
					if(wr_dist >= rd_dist){
                        for(j =0; j < (wr_dist - rd_dist); j++){
							__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
							stride +=3;
						}
						for(j =(wr_dist - rd_dist); j < wr_dist; j++){
							__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
							__builtin_prefetch (&b[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							__builtin_prefetch (&c[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							__builtin_prefetch (&d[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							__builtin_prefetch (&e[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
							stride +=3;
						}
					}else{
                        for(j =0; j < (rd_dist - wr_dist); j++){
							__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
							stride +=3;
						}
						for(j =(rd_dist - wr_dist); j < rd_dist; j++){
							__builtin_prefetch (&a[stride%(size/sizeof(double))-(rd_dist - wr_dist)], 1, 0);
							__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
							__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
							stride +=3;
						}
					}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
                        for(j =0; j < ((size/sizeof(double)) - bdist); j+=unroll){
					__builtin_prefetch (&a[stride%(size/sizeof(double))+wr_dist], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double))+rd_dist], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))+rd_dist], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))+rd_dist], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))+rd_dist], 0, 0);
						for(inner=0;inner<unroll;inner++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    stride +=3;
						}
                        }
						if(wr_dist >= rd_dist){
                        	for(j = ((size/sizeof(double)) - wr_dist); j < ((size/sizeof(double)) -rd_dist); j+=unroll){
								__builtin_prefetch (&b[stride%(size/sizeof(double))+rd_dist], 0, 0);
								__builtin_prefetch (&c[stride%(size/sizeof(double))+rd_dist], 0, 0);
								__builtin_prefetch (&d[stride%(size/sizeof(double))+rd_dist], 0, 0);
								__builtin_prefetch (&e[stride%(size/sizeof(double))+rd_dist], 0, 0);
							for(inner=0;inner<unroll;inner++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
							}
							for(j = ((size/sizeof(double)) - rd_dist); j < (size/sizeof(double)); j++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
						}else{
                        	for(j = ((size/sizeof(double)) - rd_dist); j < ((size/sizeof(double)) -wr_dist); j+=unroll){
								__builtin_prefetch (&a[stride%(size/sizeof(double))+wr_dist], 1, 0);
						for(inner=0;inner<unroll;inner++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
							}
							for(j = ((size/sizeof(double)) - rd_dist); j < (size/sizeof(double)); j++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
                 				stride +=3;
							}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			str_t += accum;
			str_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}


			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
                        for(j =0; j < rd_dist; j++){
					__builtin_prefetch (&rand_tab[j], 0, 0);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
                        for(j =0; j < ((size/sizeof(double)) - (rd_dist)); j+=unroll){
					__builtin_prefetch (&rand_tab[j+rd_dist], 0, 0);
						for(inner=0;inner<unroll;inner++)
			  a[rand_tab[j+inner]] = b[rand_tab[j+inner]] + c[rand_tab[j+inner]] + d[rand_tab[j+inner]] + e[rand_tab[j+inner]];
                        }
                        for(j = ((size/sizeof(double)) - rd_dist);j < (size/sizeof(double)); j++){
			  a[rand_tab[j]] = b[rand_tab[j]] + c[rand_tab[j]] + d[rand_tab[j]] + e[rand_tab[j]];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			rand_t += accum;
			rand_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}

			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
                        	for(j =0; j < wr_dist - rd_dist; j++){
					__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
                        	for(j =(wr_dist - rd_dist); j < wr_dist; j++){
					__builtin_prefetch (&a[stride%(size/sizeof(double))], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))-(wr_dist - rd_dist)], 0, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}else{
                        	for(j =0; j < rd_dist - wr_dist; j++){
					__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
                        	for(j =(rd_dist - wr_dist); j < rd_dist; j++){
					__builtin_prefetch (&a[stride%(size/sizeof(double))-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double))], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double))], 0, 0);
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
                        for(j =0; j < (size/sizeof(double)) - bdist; j+=unroll){
					__builtin_prefetch (&a[stride%(size/sizeof(double)) + wr_dist], 1, 0);
					__builtin_prefetch (&b[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double)) + rd_dist], 0, 0);
						for(inner=0;inner<unroll;inner++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if(((j+inner)%8 == 0)&&((j+inner) != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
                        }
				}
						if(wr_dist >= rd_dist){
                        for(j = ((size/sizeof(double)) - wr_dist); j < (size/sizeof(double)) - rd_dist; j+=unroll){
					__builtin_prefetch (&b[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&c[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&d[stride%(size/sizeof(double)) + rd_dist], 0, 0);
					__builtin_prefetch (&e[stride%(size/sizeof(double)) + rd_dist], 0, 0);
						for(inner=0;inner<unroll;inner++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if(((j+inner)%8 == 0)&&((j+inner) != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
						 for(j = (size/sizeof(double))-rd_dist; j < (size/sizeof(double)); j++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}else{
                        for(j = ((size/sizeof(double)) - rd_dist); j < (size/sizeof(double)) - wr_dist; j+=unroll){
					__builtin_prefetch (&a[stride%(size/sizeof(double)) + wr_dist], 1, 0);
						for(inner=0;inner<unroll;inner++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if(((j+inner)%8 == 0)&&((j+inner) != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
						 for(j = (size/sizeof(double))-wr_dist; j < (size/sizeof(double)); j++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
						}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			l2cache_t += accum;
                        l2cache_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}


			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist - rd_dist)){
											goto out1;
									}
									__builtin_prefetch (&a[stride], 1, 0);
							}
					}
			}
out1:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist)){
											goto out3;
									}
					__builtin_prefetch (&a[stride], 1, 0);
					__builtin_prefetch (&b[stride-(wr_dist - rd_dist)], 0, 2);
							}
					}
			}
						}else{
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist - wr_dist)){
											goto out2;
									}
					__builtin_prefetch (&b[stride], 0, 2);
						}
					}
			}
out2:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist)){
											goto out3;
									}
					__builtin_prefetch (&a[stride-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride], 0, 2);
						}
					}
			}
						}
out3:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-bdist))
											goto out5;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&a[stride+(wr_dist)], 1, 0);
					}
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&b[stride+rd_dist], 0, 2);
					}
								for(inner=0;inner<unroll;inner++)
									a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner];
							}
					}
			}
out5:
						if(wr_dist >= rd_dist){
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-rd_dist))
											goto out6;
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&b[stride+rd_dist], 0, 2);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner];
							}
					}
			}
out6:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1];
							}
					}
			}
						}else{
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-wr_dist))
											goto out7;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&a[stride+wr_dist], 1, 0);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner];
							}
					}
			}
out7:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1];
							}
					}
			}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			t_sten_t += accum;
                        t_sten_avg += ((4*size*procs*1.0E-06)/(long double)(accum - empty));
			}



			stride = 0;
			ext_flg = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist - rd_dist)){
											goto out11;
									}
									__builtin_prefetch (&a[stride], 1, 0);
							}
					}
			}
out11:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist)){
											goto out33;
									}
					__builtin_prefetch (&a[stride], 1, 0);
					__builtin_prefetch (&b[stride-(wr_dist - rd_dist)], 0, 2);
					__builtin_prefetch (&b[stride+ldim-(wr_dist - rd_dist)], 0, 2);
							}
					}
			}
						}else{
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist - wr_dist)){
											goto out22;
									}
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
						}
					}
			}
out22:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist)){
											goto out33;
									}
					__builtin_prefetch (&a[stride-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
						}
					}
			}
						}
out33:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-bdist))
											goto out55;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&a[stride+(wr_dist)], 1, 0);
					}
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&b[stride+rd_dist], 0, 2);
								__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
					}
								for(inner=0;inner<unroll;inner++)
									a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner];
							}
					}
			}
out55:
						if(wr_dist >= rd_dist){
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-rd_dist))
											goto out66;
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&b[stride+rd_dist], 0, 2);
									__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner];
							}
					}
			}
out66:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim];
							}
					}
			}
						}else{
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-wr_dist))
											goto out77;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&a[stride+wr_dist], 1, 0);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner];
							}
					}
			}
out77:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim];
							}
					}
			}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
				f_sten_t += accum;
            	f_sten_avg += ((6*size*procs*1.0E-06)/(long double)(accum - empty));
			}


			stride = 0;
			ext_flg = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist - rd_dist)){
										goto out111;
									}
											
									__builtin_prefetch (&a[stride], 1, 0);
							}
					}
			}
out111:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist)){
											goto out333;
									}
					__builtin_prefetch (&a[stride], 1, 0);
					__builtin_prefetch (&b[stride-(wr_dist - rd_dist)], 0, 2);
					__builtin_prefetch (&b[stride+ldim-(wr_dist - rd_dist)], 0, 2);
							}
					}
			}
						}else{
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist - wr_dist)){
											goto out222;
									}
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
						}
					}
			}
out222:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist)){
											goto out333;
									}
					__builtin_prefetch (&a[stride-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
						}
					}
			}
						}
out333:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-bdist))
											goto out555;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&a[stride+(wr_dist)], 1, 0);
					}
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&b[stride+rd_dist], 0, 2);
								__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
					}
								for(inner=0;inner<unroll;inner++)
									a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner] + b[stride+ldim+1+inner] + b[stride+ldim-1+inner] + b[stride-(ldim+1)+inner] + b[stride-(ldim-1)+inner];
							}
					}
			}
out555:
						if(wr_dist >= rd_dist){
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-rd_dist))
											goto out666;
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&b[stride+rd_dist], 0, 2);
									__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner] + b[stride+ldim+1+inner] + b[stride+ldim-1+inner] + b[stride-(ldim+1)+inner] + b[stride-(ldim-1)+inner];
							}
					}
			}
out666:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride+ldim+1] + b[stride+ldim-1] + b[stride-(ldim+1)] + b[stride-(ldim-1)];
							}
					}
			}
						}else{
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-wr_dist))
											goto out777;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&a[stride+wr_dist], 1, 0);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner] + b[stride+ldim+1+inner] + b[stride+ldim-1+inner] + b[stride-(ldim+1)+inner] + b[stride-(ldim-1)+inner];
							}
					}
			}
out777:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride+ldim+1] + b[stride+ldim-1] + b[stride-(ldim+1)] + b[stride-(ldim-1)];
							}
					}
			}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
				n_sten_t += accum;
            	n_sten_avg += ((10*size*procs*1.0E-06)/(long double)(accum - empty));
			}


			stride = 0;
			ext_flg = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist - rd_dist))
									{
											goto out1111;
									}
									__builtin_prefetch (&a[stride], 1, 0);
							}
					}
			}
out1111:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist))
									{
											goto out3333;
									}
					__builtin_prefetch (&a[stride], 1, 0);
					__builtin_prefetch (&b[stride-(wr_dist - rd_dist)], 0, 2);
					__builtin_prefetch (&b[stride+ldim-(wr_dist - rd_dist)], 0, 2);
					__builtin_prefetch (&b[stride+(ldim*ldim)-(wr_dist - rd_dist)], 0, 2);
							}
					}
			}
						}else{
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist - wr_dist))
									{
											goto out2222;
									}
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
					__builtin_prefetch (&b[stride+(ldim*ldim)], 0, 2);
						}
					}
			}
out2222:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist))
									{
											goto out3333;
									}
					__builtin_prefetch (&a[stride-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
					__builtin_prefetch (&b[stride+(ldim*ldim)], 0, 2);
						}
					}
			}
						}
out3333:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-bdist))
											goto out5555;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&a[stride+(wr_dist)], 1, 0);
					}
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&b[stride+rd_dist], 0, 2);
								__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
								__builtin_prefetch (&b[stride+(ldim*ldim)+rd_dist], 0, 2);
					}
								for(inner=0;inner<unroll;inner++)
									a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner] + b[stride + ldim*ldim+inner] + b[stride - ldim*ldim+inner];
							}
					}
			}
out5555:
						if(wr_dist >= rd_dist){
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-rd_dist))
											goto out6666;
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&b[stride+rd_dist], 0, 2);
									__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
									__builtin_prefetch (&b[stride+(ldim*ldim)+rd_dist], 0, 2);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner] + b[stride + ldim*ldim+inner] + b[stride - ldim*ldim+inner];
							}
					}
			}
out6666:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride + ldim*ldim] + b[stride - ldim*ldim];
							}
					}
			}
						}else{
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-wr_dist))
											goto out7777;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&a[stride+wr_dist], 1, 0);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+1+inner] + b[stride-1+inner] + b[stride+ldim+inner] + b[stride-ldim+inner] + b[stride + ldim*ldim+inner] + b[stride - ldim*ldim+inner];
							}
					}
			}
out7777:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride + ldim*ldim] + b[stride - ldim*ldim];
							}
					}
			}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
				s_sten_t += accum;
            	s_sten_avg += ((8*size*procs*1.0E-06)/(long double)(accum - empty));
			}


			stride = 0;
			ext_flg = 0;
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &obegin);
//#pragma omp parallel for
						if(wr_dist >= rd_dist){
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist - rd_dist))
									{
											goto out11111;
									}
									__builtin_prefetch (&a[stride], 1, 0);
							}
					}
			}
out11111:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (wr_dist))
									{
										goto out33333;
									}
					__builtin_prefetch (&a[stride], 1, 0);
					__builtin_prefetch (&b[stride-(wr_dist - rd_dist)], 0, 2);
					__builtin_prefetch (&b[stride+ldim-(wr_dist - rd_dist)], 0, 2);
					__builtin_prefetch (&b[stride+(ldim*ldim)-(wr_dist - rd_dist)], 0, 2);
							}
					}
			}
						}else{
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist - wr_dist))
									{
											goto out22222;
									}
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
					__builtin_prefetch (&b[stride+(ldim*ldim)], 0, 2);
						}
					}
			}
out22222:
			ext_flg = 0;
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									if(stride > (rd_dist))
									{
											goto out33333;
									}
					__builtin_prefetch (&a[stride-(rd_dist - wr_dist)], 1, 0);
					__builtin_prefetch (&b[stride], 0, 2);
					__builtin_prefetch (&b[stride+ldim], 0, 2);
					__builtin_prefetch (&b[stride+(ldim*ldim)], 0, 2);
						}
					}
			}
						}
out33333:
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &ostop);
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ( (ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-bdist))
											goto out55555;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&a[stride+(wr_dist)], 1, 0);
					}
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
								__builtin_prefetch (&b[stride+rd_dist], 0, 2);
								__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
								__builtin_prefetch (&b[stride+(ldim*ldim)+rd_dist], 0, 2);
					}
								for(inner=0;inner<unroll;inner++)
									a[stride+inner] = b[stride+inner] + b[stride+inner+1] + b[stride+inner-1] + b[stride+inner+ldim] + b[stride+inner-ldim] + b[stride+inner+ldim+1] + b[stride+inner+ldim-1] + b[stride+inner-(ldim+1)] + b[stride+inner-(ldim-1)] + b[stride+inner+(ldim*ldim)] + b[stride+inner-(ldim*ldim)] + b[stride+inner+(ldim*ldim)+1] + b[stride+inner-(ldim*ldim)+1] + b[stride+inner+(ldim*ldim)-1] + b[stride+inner-(ldim*ldim)-1] + b[stride+inner+(ldim*ldim)+ldim] + b[stride+inner-(ldim*ldim)+ldim] + b[stride+inner+(ldim*ldim)-ldim] + b[stride+inner-(ldim*ldim)-ldim] + b[stride+inner+(ldim*ldim)+ldim+1] + b[stride+inner-(ldim*ldim)+ldim+1] + b[stride+inner+(ldim*ldim)-ldim+1] + b[stride+inner-(ldim*ldim)-ldim+1] + b[stride+inner+(ldim*ldim)+ldim-1] + b[stride+inner-(ldim*ldim)+ldim-1] + b[stride+inner+(ldim*ldim)-ldim-1] + b[stride+inner-(ldim*ldim)-ldim-1];
							}
					}
			}
out55555:
						if(wr_dist >= rd_dist){
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-rd_dist))
											goto out66666;
					if((rd_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&b[stride+rd_dist], 0, 2);
									__builtin_prefetch (&b[stride+ldim+rd_dist], 0, 2);
									__builtin_prefetch (&b[stride+(ldim*ldim)+rd_dist], 0, 2);
					}
									for(inner=0;inner<unroll;inner++)
										a[stride+inner] = b[stride+inner] + b[stride+inner+1] + b[stride+inner-1] + b[stride+inner+ldim] + b[stride+inner-ldim] + b[stride+inner+ldim+1] + b[stride+inner+ldim-1] + b[stride+inner-(ldim+1)] + b[stride+inner-(ldim-1)] + b[stride+inner+(ldim*ldim)] + b[stride+inner-(ldim*ldim)] + b[stride+inner+(ldim*ldim)+1] + b[stride+inner-(ldim*ldim)+1] + b[stride+inner+(ldim*ldim)-1] + b[stride+inner-(ldim*ldim)-1] + b[stride+inner+(ldim*ldim)+ldim] + b[stride+inner-(ldim*ldim)+ldim] + b[stride+inner+(ldim*ldim)-ldim] + b[stride+inner-(ldim*ldim)-ldim] + b[stride+inner+(ldim*ldim)+ldim+1] + b[stride+inner-(ldim*ldim)+ldim+1] + b[stride+inner+(ldim*ldim)-ldim+1] + b[stride+inner-(ldim*ldim)-ldim+1] + b[stride+inner+(ldim*ldim)+ldim-1] + b[stride+inner-(ldim*ldim)+ldim-1] + b[stride+inner+(ldim*ldim)-ldim-1] + b[stride+inner-(ldim*ldim)-ldim-1];
							}
					}
			}
out66666:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride+ldim+1] + b[stride+ldim-1] + b[stride-(ldim+1)] + b[stride-(ldim-1)] + b[stride+(ldim*ldim)] + b[stride-(ldim*ldim)] + b[stride+(ldim*ldim)+1] + b[stride-(ldim*ldim)+1] + b[stride+(ldim*ldim)-1] + b[stride-(ldim*ldim)-1] + b[stride+(ldim*ldim)+ldim] + b[stride-(ldim*ldim)+ldim] + b[stride+(ldim*ldim)-ldim] + b[stride-(ldim*ldim)-ldim] + b[stride+(ldim*ldim)+ldim+1] + b[stride-(ldim*ldim)+ldim+1] + b[stride+(ldim*ldim)-ldim+1] + b[stride-(ldim*ldim)-ldim+1] + b[stride+(ldim*ldim)+ldim-1] + b[stride-(ldim*ldim)+ldim-1] + b[stride+(ldim*ldim)-ldim-1] + b[stride-(ldim*ldim)-ldim-1];
							}
					}
			}
						}else{
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k+=unroll){
									stride = l + j + k;
									if(stride > ((ldim + 1)*(ldim + 1)*(ldim)-wr_dist))
											goto out77777;
					if((wr_dist < pg_sz)||(stride%(pg_sz/sizeof(double)) == 0)){
									__builtin_prefetch (&a[stride+wr_dist], 1, 0);
					}
									for(inner=0;inner<unroll;inner++)
									a[stride+inner] = b[stride+inner] + b[stride+inner+1] + b[stride+inner-1] + b[stride+inner+ldim] + b[stride+inner-ldim] + b[stride+inner+ldim+1] + b[stride+inner+ldim-1] + b[stride+inner-(ldim+1)] + b[stride+inner-(ldim-1)] + b[stride+inner+(ldim*ldim)] + b[stride+inner-(ldim*ldim)] + b[stride+inner+(ldim*ldim)+1] + b[stride+inner-(ldim*ldim)+1] + b[stride+inner+(ldim*ldim)-1] + b[stride+inner-(ldim*ldim)-1] + b[stride+inner+(ldim*ldim)+ldim] + b[stride+inner-(ldim*ldim)+ldim] + b[stride+inner+(ldim*ldim)-ldim] + b[stride+inner-(ldim*ldim)-ldim] + b[stride+inner+(ldim*ldim)+ldim+1] + b[stride+inner-(ldim*ldim)+ldim+1] + b[stride+inner+(ldim*ldim)-ldim+1] + b[stride+inner-(ldim*ldim)-ldim+1] + b[stride+inner+(ldim*ldim)+ldim-1] + b[stride+inner-(ldim*ldim)+ldim-1] + b[stride+inner+(ldim*ldim)-ldim-1] + b[stride+inner-(ldim*ldim)-ldim-1];
							}
					}
			}
out77777:
			for(l = l; l < ((ldim + 1)*(ldim + 1)*(ldim)); l+=(ldim + 1)*(ldim + 1)){
					for(j = j; j < ((ldim +1)*(ldim + 1) - (ldim + 1)); j += (ldim + 1)){
							for(k = k; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride+ldim+1] + b[stride+ldim-1] + b[stride-(ldim+1)] + b[stride-(ldim-1)] + b[stride+(ldim*ldim)] + b[stride-(ldim*ldim)] + b[stride+(ldim*ldim)+1] + b[stride-(ldim*ldim)+1] + b[stride+(ldim*ldim)-1] + b[stride-(ldim*ldim)-1] + b[stride+(ldim*ldim)+ldim] + b[stride-(ldim*ldim)+ldim] + b[stride+(ldim*ldim)-ldim] + b[stride-(ldim*ldim)-ldim] + b[stride+(ldim*ldim)+ldim+1] + b[stride-(ldim*ldim)+ldim+1] + b[stride+(ldim*ldim)-ldim+1] + b[stride-(ldim*ldim)-ldim+1] + b[stride+(ldim*ldim)+ldim-1] + b[stride-(ldim*ldim)+ldim-1] + b[stride+(ldim*ldim)-ldim-1] + b[stride-(ldim*ldim)-ldim-1];
							}
					}
			}
						}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
				t7_sten_t += accum;
            	t7_sten_avg += ((28*size*procs*1.0E-06)/(long double)(accum - empty));
			}
*/

			numa_free(a, size);
			numa_free(b, size);
			numa_free(c, size);
			numa_free(d, size);
			numa_free(e, size);
		}
		if(rank == 0)
		{
			printf("%d %d-%d %d %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf\n", procs, wr_dist*sizeof(double), rd_dist*sizeof(double), unroll, wr_only_avg/10, owfr_avg/10, str_avg/10, rand_avg/10, l2cache_avg/10, t_sten_avg/10, f_sten_avg/10, s_sten_avg/10, n_sten_avg/10, t7_sten_avg/10);
		}
				unroll*=2;
				}
			rd_dist*=2;
		}
		if(rank == 0)
		printf("\n");
	wr_dist*=2;
}
/*if(rank == 0){
		struct numa_node_bw * node_bw = (struct numa_node_bw *)malloc(sizeof(struct numa_node_bw));
		node_bw->numa_id = numa_node_ids[i];
		node_bw->wr_only_avg = wr_only_avg/10;
		node_bw->owfr_avg = owfr_avg/10;
		node_bw->str_avg = str_avg/10;
		node_bw->rand_avg = rand_avg/10;
		node_bw->l2cache_avg = l2cache_avg/10;
		node_bw->t_sten_avg = t_sten_avg/10;
		node_bw->f_sten_avg = f_sten_avg/10;
		node_bw->n_sten_avg = n_sten_avg/10;
		node_bw->s_sten_avg = s_sten_avg/10;
		node_bw->t7_sten_avg = t7_sten_avg/10;
		node_bw->wr_only_t = wr_only_t/10;
		node_bw->l2cache_t = l2cache_t/10;
		node_bw->rand_t = rand_t/10;
		node_bw->str_t = str_t/10;
		node_bw->owfr_t = owfr_t/10;
		node_bw->t_sten_t = t_sten_t/10;
		node_bw->f_sten_t = f_sten_t/10;
		node_bw->n_sten_t = n_sten_t/10;
		node_bw->s_sten_t = s_sten_t/10;
		node_bw->t7_sten_t = t7_sten_t/10;
		node_bw->next = NULL;
		}*/
		i+=2;
		if(rank == 0)
		{
				printf("\n\n");
		}
	}
//	if(rank == 0){
//	write_config_file();
//	}
	free(rand_tab);
}
