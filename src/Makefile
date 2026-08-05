CC=g++
CFLAGS= -I. -std=c++17 -O3 -march=native -flto=auto -fopenmp -fno-math-errno
# Perhaps you need to add a path to omp headers like -I/path/to/omp/headers
# OSX clang does not support omp by default, you need to install library and compiler separately (e.g. homebrew)
# For OSX you need to change -march to -mcpu or remove it
LDFLAGS= -lz
Sanity: calc_true_variation_parallel_prior_mu_sigma.o ReadInputFiles.o FitFrac.o FukushimaLambertW.o Digamma_Trigamma.o npy_writer.o
	$(CC) $(CFLAGS) -o ../bin/Sanity calc_true_variation_parallel_prior_mu_sigma.o ReadInputFiles.o FitFrac.o FukushimaLambertW.o Digamma_Trigamma.o npy_writer.o $(LDFLAGS)

calc_true_variation_parallel_prior_mu_sigma.o: calc_true_variation_parallel_prior_mu_sigma.cpp ReadInputFiles.h FitFrac.h FukushimaLambertW.h Digamma_Trigamma.h Writer.hpp npy_writer.hpp
	$(CC) $(CFLAGS) -c calc_true_variation_parallel_prior_mu_sigma.cpp

npy_writer.o: npy_writer.cpp npy_writer.hpp
	$(CC) $(CFLAGS) -c npy_writer.cpp

ReadInputFiles.o: ReadInputFiles.cpp ReadInputFiles.h
	$(CC) $(CFLAGS) -c ReadInputFiles.cpp

FitFrac.o: FitFrac.cpp FitFrac.h FukushimaLambertW.h
	$(CC) $(CFLAGS) -c FitFrac.cpp

FukushimaLambertW.o: FukushimaLambertW.cc FukushimaLambertW.h
	$(CC) $(CFLAGS) -c FukushimaLambertW.cc

Digamma_Trigamma.o: Digamma_Trigamma.cpp Digamma_Trigamma.h
	$(CC) $(CFLAGS) -c Digamma_Trigamma.cpp

Sanity_distance: compute_distance.o
	$(CC) $(CFLAGS) -o ../bin/Sanity_distance compute_distance.o

compute_distance.o: compute_distance.cpp
	$(CC) $(CFLAGS) -c compute_distance.cpp

Sanity_gene_correlation:  compute_gene_correlations.o
	$(CC) $(CFLAGS) -o ../bin/Sanity_gene_correlation compute_gene_correlations.o

compute_gene_correlation.o: compute_gene_correlations.cpp
	$(CC) $(CFLAGS) -c compute_gene_correlations.cpp

clean:
	rm -f *.o

