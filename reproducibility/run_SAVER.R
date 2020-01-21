# load libraries
library(SAVER)

# Define input file and output directory
dataset <- "path/to/UMI_count_table.txt"
outdir <- "path/to/output/directory"
out.path<- paste(outdir, 'SAVER_normalization.txt', sep='/')

# Read count matrix
raw.data <- read.table(dataset, header = TRUE, row.names = 1)
count <- as.matrix(raw.data)

# Run SAVER
data.saver <- saver(count, ncores = 8, estimates.only = TRUE)

# Write output matrix
write.table(data.saver, file=out.path, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
