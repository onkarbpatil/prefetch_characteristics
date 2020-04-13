CC = mpicc
LDFLAGS = -O1 -fprefetch-loop-arrays -DLIKWID_PERFMON -fopenmp -lnuma -lm -L/usr/local/lib -I/usr/local/include -llikwid

numatest_mpi: test3.c numatest_mpi_np.c numatest_mpi.h
		$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
	clean:
		rm -f $(obj) numatest_mpi
