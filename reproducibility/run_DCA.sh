#!/bin/bash

# Activate DCA conda environment
source ~/miniconda3/bin/activate DCA

# Define input file and output directory
DATASET=path/to/UMI_count_table.txt
OUTDIR=path/to/output/directory

# Run DCA
dca $DATASET $OUTDIR --threads 8
