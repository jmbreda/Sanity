# Import libraries
import magic
import pandas as pd
import numpy as np
import os

# Define input file and output directory
dataset = "path/to/UMI_count_table.txt"
outdir = "path/to/output/directory"

# Read count matrix
X = pd.read_csv(dataset,sep='\t',index_col=0)

# Remove genes summing to 0 counts
X = X.loc[X.sum(axis=1) > 0,:]

# Transpose to cell X gene
X = X.T

# Run MAGIC
magic_op = magic.MAGIC()
magic_op.set_params(n_pca=min(X.shape[0],100))
X_magic = magic_op.fit_transform(X,genes='all_genes')

# Transpose back to gene X cell
X_magic = X_magic.T

# Write output matrix
X_magic.to_csv(outdir + '/magic_normalization.txt',sep='\t')
