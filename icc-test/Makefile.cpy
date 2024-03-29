CC = mpicc
LDFLAGS = -O3 -qopt-prefetch -lnuma -lm

numatest_mpi: test3.c numatest_mpi_np.c numatest_mpi.h
		$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
	clean:
		rm -f $(obj) numatest_mpi
