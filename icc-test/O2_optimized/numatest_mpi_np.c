#include "numatest_mpi.h"

int mem_types;
int max_node;
int numt;
int total_numa_nodes = 0;
int * numa_node_ids;
struct bitmask * numa_nodes;
char ** mem_tech;
long double * means;
int * cluster_sizes;



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
	int rs = 0;
	int z = 0;
	int dist = 0;
	int rd_dist;
			int wr_dist;
	int * rand_tab;
	rand_tab = (int*)malloc(mbs*sizeof(int));
	double *a, *b, *c, *d, *e, *f, *g, *h;
	double **aa, **bb, **cc;
	clock_t start, end;
	struct timespec begin, stop;
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

//#ifdef _OPENMP
//#pragma omp parallel private(numt)
    {
//    numt = omp_get_num_threads();
    }
//#endif
  	i = 0;
	while(i < total_numa_nodes){
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
		for( iters = 0; iters < 10; iters++)
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
			//LIKWID_MARKER_START("wr-only");
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j = 0;j < (size/sizeof(double));j++){
				a[j] = 1.0;
				b[j] = 2.0;
				c[j] = 3.0;
				d[j] = 4.0;
				e[j] = 5.0;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("wr-only");
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			wr_only_t += accum;
			wr_only_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}

			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("1w4r");
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
			for(j =0; j < (size/sizeof(double)); j++){
                            a[j] = c[j] + d[j] + e[j] + b[j];
            }
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("1w4r");
			if(rank == 0){
			accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			owfr_t += accum;
			owfr_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
/*
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("str");
			clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
								a[stride%(size/sizeof(double))] = c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + b[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    stride +=3;
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("str");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			str_t += accum;
			str_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}

			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("rand");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
			  a[rand_tab[j]] = b[rand_tab[j]] + c[rand_tab[j]] + d[rand_tab[j]] + e[rand_tab[j]];
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("rand");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			rand_t += accum;
			rand_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}

			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("l3");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
//#pragma omp parallel for
                        for(j =0; j < (size/sizeof(double)); j++){
                            a[stride%(size/sizeof(double))] = b[stride%(size/sizeof(double))] + c[stride%(size/sizeof(double))] + d[stride%(size/sizeof(double))] + e[stride%(size/sizeof(double))];
			    if((j%8 == 0)&&(j != 0))
				stride = j*4757914; //65536 for KNL
			    else
				stride++;
                        }
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("l3");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			l2cache_t += accum;
                        l2cache_avg += ((5*size*procs*1.0E-06)/(long double)(accum - empty));
			}
*/
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("3sten");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ((ldim +1)*(ldim + 1) - (ldim+1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1];
							}
					}
			}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("3sten");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			t_sten_t += accum;
			t_sten_avg += ((4*size*procs*1.0E-06)/(long double)(accum - empty));
			}
			
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("5sten");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ((ldim +1)*(ldim + 1) - (ldim+1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim];
							}
					}
			}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("5sten");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			f_sten_t += accum;
			f_sten_avg += ((6*size*procs*1.0E-06)/(long double)(accum - empty));
			}
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("7sten");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = ldim*ldim; l < ((ldim+1)*(ldim+1)*ldim);l+=(ldim+1)*(ldim+1)){
					for(j = ldim; j < (((ldim+1)*(ldim+1)) - ldim); j += (ldim+1)){
							for(k = 1; k < ldim; k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride + ldim*ldim] + b[stride - ldim*ldim];
							}
					}
			}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("7sten");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			s_sten_t += accum;
			s_sten_avg += ((8*size*procs*1.0E-06)/(long double)(accum - empty));
			}
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("9sten");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ((ldim +1)*(ldim + 1) - (ldim+1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride+ldim+1] + b[stride+ldim-1] + b[stride-(ldim+1)] + b[stride-(ldim-1)];
							}
					}
			}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("9sten");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			n_sten_t += accum;
			n_sten_avg += ((10*size*procs*1.0E-06)/(long double)(accum - empty));
			}
			stride = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			//LIKWID_MARKER_START("27sten");
                        clock_gettime( CLOCK_MONOTONIC, &begin);
			for(l = (ldim + 1)*(ldim + 1); l < ((ldim + 1)*(ldim + 1)*(ldim));l+=(ldim + 1)*(ldim + 1)){
					for(j = (ldim+1); j < ((ldim +1)*(ldim + 1) - (ldim+1)); j += (ldim + 1)){
							for(k = 1; k < (ldim); k++){
									stride = l + j + k;
									a[stride] = b[stride] + b[stride+1] + b[stride-1] + b[stride+ldim] + b[stride-ldim] + b[stride+ldim+1] + b[stride+ldim-1] + b[stride-(ldim+1)] + b[stride-(ldim-1)] + b[stride+(ldim*ldim)] + b[stride-(ldim*ldim)] + b[stride+(ldim*ldim)+1] + b[stride-(ldim*ldim)+1] + b[stride+(ldim*ldim)-1] + b[stride-(ldim*ldim)-1] + b[stride+(ldim*ldim)+ldim] + b[stride-(ldim*ldim)+ldim] + b[stride+(ldim*ldim)-ldim] + b[stride-(ldim*ldim)-ldim] + b[stride+(ldim*ldim)+ldim+1] + b[stride-(ldim*ldim)+ldim+1] + b[stride+(ldim*ldim)-ldim+1] + b[stride-(ldim*ldim)-ldim+1] + b[stride+(ldim*ldim)+ldim-1] + b[stride-(ldim*ldim)+ldim-1] + b[stride+(ldim*ldim)-ldim-1] + b[stride-(ldim*ldim)-ldim-1];
							}
					}
			}
			MPI_Barrier(MPI_COMM_WORLD);
                        clock_gettime( CLOCK_MONOTONIC, &stop);
			//LIKWID_MARKER_STOP("27sten");
			if(rank == 0){
                        accum = ( stop.tv_sec - begin.tv_sec ) + (long double)( stop.tv_nsec - begin.tv_nsec ) / (long double)BILLION;
			t7_sten_t += accum;
			t7_sten_avg += ((28*size*procs*1.0E-06)/(long double)(accum - empty));
			}
			
			numa_free(a, size);
			numa_free(b, size);
			numa_free(c, size);
			numa_free(d, size);
			numa_free(e, size);
		}
/*		if(rank == 0){
		struct numa_node_bw * node_bw = (struct numa_node_bw *)malloc(sizeof(struct numa_node_bw));
		node_bw->numa_id = numa_node_ids[i];
		node_bw->wr_only_avg = wr_only_avg;
		node_bw->owfr_avg = owfr_avg;
		node_bw->str_avg = str_avg;
		node_bw->rand_avg = rand_avg;
		node_bw->l2cache_avg = l2cache_avg;
		node_bw->t_sten_avg = t_sten_avg;
		node_bw->f_sten_avg = f_sten_avg;
		node_bw->n_sten_avg = n_sten_avg;
		node_bw->s_sten_avg = s_sten_avg;
		node_bw->t7_sten_avg = t7_sten_avg;
		node_bw->wr_only_t = wr_only_t;
		node_bw->l2cache_t = l2cache_t;
		node_bw->rand_t = rand_t;
		node_bw->str_t = str_t;
		node_bw->owfr_t = owfr_t;
		node_bw->t_sten_t = t_sten_t;
		node_bw->f_sten_t = f_sten_t;
		node_bw->n_sten_t = n_sten_t;
		node_bw->s_sten_t = s_sten_t;
		node_bw->t7_sten_t = t7_sten_t;
		node_bw->next = NULL;
		if(numa_node_list == NULL){
			numa_node_list = node_bw;
			numa_list_head = numa_node_list;
		}
		else{
			sort_list(node_bw);
		}*/
		if(rank == 0)
		{
			printf("%d %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf\n\n", i, wr_only_avg, owfr_avg, str_avg, rand_avg, l2cache_avg, t_sten_avg, f_sten_avg, s_sten_avg, n_sten_avg, t7_sten_avg);
		}
		i+=2;
	}
	free(rand_tab);
}
