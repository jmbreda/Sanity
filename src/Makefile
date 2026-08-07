UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
    # Apple's clang (installed as g++/c++) does not support -fopenmp, so we
    # need a real GCC. Install one with `brew install gcc`, or pass your own
    # compiler on the command line, e.g. `make CC=g++-mp-13` for MacPorts.
    ifneq ($(origin CC),command line)
        GCC_PREFIX := $(shell brew --prefix gcc 2>/dev/null)
        ifeq ($(GCC_PREFIX),)
            $(error GCC not found. Install it with: brew install gcc)
        endif
        CC := $(lastword $(sort $(wildcard $(GCC_PREFIX)/bin/g++-*)))
        ifeq ($(CC),)
            $(error No g++-* binary found under $(GCC_PREFIX)/bin. Install GCC with: brew install gcc)
        endif
    endif
    ifeq ($(UNAME_M),arm64)
        ARCH_FLAG := -mcpu=native
    else
        ARCH_FLAG := -march=native
    endif
else
    CC=g++
    ARCH_FLAG := -march=native
endif

CFLAGS= -I. -std=c++17 -O3 $(ARCH_FLAG) -flto=auto -fopenmp -fno-math-errno
LDFLAGS= -lz
Sanity: calc_true_variation_parallel_prior_mu_sigma.o ReadInputFiles.o FitFrac.o FukushimaLambertW.o Digamma_Trigamma.o npy_writer.o
	$(CC) $(CFLAGS) -o ../bin/Sanity calc_true_variation_parallel_prior_mu_sigma.o ReadInputFiles.o FitFrac.o FukushimaLambertW.o Digamma_Trigamma.o npy_writer.o $(LDFLAGS)

calc_true_variation_parallel_prior_mu_sigma.o: calc_true_variation_parallel_prior_mu_sigma.cpp ReadInputFiles.h FitFrac.h FukushimaLambertW.h Digamma_Trigamma.h Writer.hpp npy_writer.hpp Version.h
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

compute_distance.o: compute_distance.cpp Version.h
	$(CC) $(CFLAGS) -c compute_distance.cpp

Sanity_gene_correlation:  compute_gene_correlations.o
	$(CC) $(CFLAGS) -o ../bin/Sanity_gene_correlation compute_gene_correlations.o

compute_gene_correlations.o: compute_gene_correlations.cpp Version.h
	$(CC) $(CFLAGS) -c compute_gene_correlations.cpp

.PHONY: clean

clean:
	rm -f *.o
	rm -f ../bin/Sanity ../bin/Sanity_distance ../bin/Sanity_gene_correlation

