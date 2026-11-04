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

SANITY_OBJS := sanity_main.o sanity_options.o gene_fit.o fit_frac.o \
    read_input_files.o dense_row_reader.o npy_writer.o FukushimaLambertW.o digamma_trigamma.o

Sanity: $(SANITY_OBJS) | $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/Sanity $(SANITY_OBJS) $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

sanity_main.o: sanity_main.cpp gene_fit.hpp logging.hpp sanity_options.hpp Version.h $(IO)/format_output.hpp $(IO)/read_input_files.hpp $(IO)/dense_row_reader.hpp $(IO)/writer.hpp $(IO)/npy_writer.hpp
	$(CC) $(CFLAGS) -c sanity_main.cpp

sanity_options.o: sanity_options.cpp sanity_options.hpp gene_fit.hpp logging.hpp
	$(CC) $(CFLAGS) -c sanity_options.cpp

gene_fit.o: gene_fit.cpp gene_fit.hpp fit_frac.hpp $(MATH)/digamma_trigamma.hpp
	$(CC) $(CFLAGS) -c gene_fit.cpp

fit_frac.o: fit_frac.cpp fit_frac.hpp $(MATH)/FukushimaLambertW.h
	$(CC) $(CFLAGS) -c fit_frac.cpp

dense_row_reader.o: $(IO)/dense_row_reader.cpp $(IO)/dense_row_reader.hpp $(IO)/read_input_files.hpp
	$(CC) $(CFLAGS) -c $(IO)/dense_row_reader.cpp

npy_writer.o: $(IO)/npy_writer.cpp $(IO)/npy_writer.hpp
	$(CC) $(CFLAGS) -c $(IO)/npy_writer.cpp

read_input_files.o: $(IO)/read_input_files.cpp $(IO)/read_input_files.hpp $(IO)/fast_mtx_reader.hpp
	$(CC) $(CFLAGS) -c $(IO)/read_input_files.cpp

FukushimaLambertW.o: $(MATH)/FukushimaLambertW.cc $(MATH)/FukushimaLambertW.h
	$(CC) $(CFLAGS) -c $(MATH)/FukushimaLambertW.cc

digamma_trigamma.o: $(MATH)/digamma_trigamma.cpp $(MATH)/digamma_trigamma.hpp
	$(CC) $(CFLAGS) -c $(MATH)/digamma_trigamma.cpp

Sanity_distance: sanity_distance.o | $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/Sanity_distance sanity_distance.o

sanity_distance.o: sanity_distance.cpp Version.h
	$(CC) $(CFLAGS) -c sanity_distance.cpp

Sanity_gene_correlation:  sanity_gene_correlation.o | $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/Sanity_gene_correlation sanity_gene_correlation.o

sanity_gene_correlation.o: sanity_gene_correlation.cpp Version.h
	$(CC) $(CFLAGS) -c sanity_gene_correlation.cpp

.PHONY: clean

clean:
	rm -f *.o
	rm -f $(BINDIR)/Sanity $(BINDIR)/Sanity_distance $(BINDIR)/Sanity_gene_correlation
