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
BINDIR := ../bin

IO := helpers/io
MATH := helpers/math

SANITY_OBJS := calc_true_variation_parallel_prior_mu_sigma.o FitFrac.o \
    ReadInputFiles.o DenseRowReader.o npy_writer.o FukushimaLambertW.o Digamma_Trigamma.o

Sanity: $(SANITY_OBJS) | $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/Sanity $(SANITY_OBJS) $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

calc_true_variation_parallel_prior_mu_sigma.o: calc_true_variation_parallel_prior_mu_sigma.cpp FitFrac.h Version.h $(IO)/FormatOutput.hpp $(IO)/ReadInputFiles.h $(IO)/DenseRowReader.hpp $(IO)/Writer.hpp $(IO)/npy_writer.hpp $(MATH)/Digamma_Trigamma.h
	$(CC) $(CFLAGS) -c calc_true_variation_parallel_prior_mu_sigma.cpp

FitFrac.o: FitFrac.cpp FitFrac.h $(MATH)/FukushimaLambertW.h
	$(CC) $(CFLAGS) -c FitFrac.cpp

DenseRowReader.o: $(IO)/DenseRowReader.cpp $(IO)/DenseRowReader.hpp $(IO)/ReadInputFiles.h
	$(CC) $(CFLAGS) -c $(IO)/DenseRowReader.cpp

npy_writer.o: $(IO)/npy_writer.cpp $(IO)/npy_writer.hpp
	$(CC) $(CFLAGS) -c $(IO)/npy_writer.cpp

ReadInputFiles.o: $(IO)/ReadInputFiles.cpp $(IO)/ReadInputFiles.h $(IO)/FastMtxReader.hpp
	$(CC) $(CFLAGS) -c $(IO)/ReadInputFiles.cpp

FukushimaLambertW.o: $(MATH)/FukushimaLambertW.cc $(MATH)/FukushimaLambertW.h
	$(CC) $(CFLAGS) -c $(MATH)/FukushimaLambertW.cc

Digamma_Trigamma.o: $(MATH)/Digamma_Trigamma.cpp $(MATH)/Digamma_Trigamma.h
	$(CC) $(CFLAGS) -c $(MATH)/Digamma_Trigamma.cpp

Sanity_distance: compute_distance.o | $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/Sanity_distance compute_distance.o

compute_distance.o: compute_distance.cpp Version.h
	$(CC) $(CFLAGS) -c compute_distance.cpp

Sanity_gene_correlation:  compute_gene_correlations.o | $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/Sanity_gene_correlation compute_gene_correlations.o

compute_gene_correlations.o: compute_gene_correlations.cpp Version.h
	$(CC) $(CFLAGS) -c compute_gene_correlations.cpp

.PHONY: clean

clean:
	rm -f *.o
	rm -f $(BINDIR)/Sanity $(BINDIR)/Sanity_distance $(BINDIR)/Sanity_gene_correlation
