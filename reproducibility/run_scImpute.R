# load libraries
library(scImpute)

# Define input file and output directory
dataset <- "path/to/UMI_count_table.txt"
outdir <- "path/to/output/directory/"

# run scImpute
scimpute(count_path = dataset, infile = "txt", outfile = "txt", out_dir=outdir, labeled = FALSE, drop_thre = 0.5, Kcluster = 1, ncore = 8)
