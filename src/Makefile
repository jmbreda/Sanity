CC=g++
CFLAGS= -I. -O2 -fopenmp -std=c++11

TrueVar: calc_true_variation_parallel_prior_mu_sigma.o FitFrac.o FukushimaLambertW.o FitVarPrior_NewtonRaphson_Backtracking.o Digamma_Trigamma.o
	$(CC) $(CFLAGS) -o ../bin/TrueVar calc_true_variation_parallel_prior_mu_sigma.o FitFrac.o FukushimaLambertW.o FitVarPrior_NewtonRaphson_Backtracking.o Digamma_Trigamma.o

calc_true_variation_parallel_prior_mu_sigma.o: calc_true_variation_parallel_prior_mu_sigma.cpp FitFrac.h FukushimaLambertW.h FitVarPrior_NewtonRaphson_Backtracking.h Digamma_Trigamma.h
	$(CC) $(CFLAGS) -c calc_true_variation_parallel_prior_mu_sigma.cpp

FitFrac.o: FitFrac.cpp FitFrac.h FukushimaLambertW.h
	$(CC) $(CFLAGS) -c FitFrac.cpp

FukushimaLambertW.o: FukushimaLambertW.cc FukushimaLambertW.h
	$(CC) $(CFLAGS) -c FukushimaLambertW.cc

FitVarPrior_NewtonRaphson_Backtracking.o: FitVarPrior_NewtonRaphson_Backtracking.cpp FitVarPrior_NewtonRaphson_Backtracking.h
	$(CC) $(CFLAGS) -c FitVarPrior_NewtonRaphson_Backtracking.cpp

Digamma_Trigamma.o: Digamma_Trigamma.cpp Digamma_Trigamma.h
	$(CC) $(CFLAGS) -c Digamma_Trigamma.cpp

clean:
	rm -f *.o

