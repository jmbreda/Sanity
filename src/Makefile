CC=g++
CFLAGS= -I. -O2 -ffast-math -O3 -fopenmp -std=c++11

Sanity: calc_true_variation_parallel_prior_mu_sigma.o ReadInputFiles.o FitFrac.o FukushimaLambertW.o Digamma_Trigamma.o
	$(CC) $(CFLAGS) -o ../bin/Sanity calc_true_variation_parallel_prior_mu_sigma.o ReadInputFiles.o FitFrac.o FukushimaLambertW.o Digamma_Trigamma.o

calc_true_variation_parallel_prior_mu_sigma.o: calc_true_variation_parallel_prior_mu_sigma.cpp ReadInputFiles.h FitFrac.h FukushimaLambertW.h Digamma_Trigamma.h
	$(CC) $(CFLAGS) -c calc_true_variation_parallel_prior_mu_sigma.cpp

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

