CC = mpicc
LDFLAGS = -fprefetch-loop-arrays -fopenmp -lnuma -lm

numatest_mpi: test3.c numatest_mpi_np.c numatest_mpi.h
		$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
	clean:
		rm -f $(obj) numatest_mpi
