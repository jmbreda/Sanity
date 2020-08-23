# Define input file and output directory
dataset <- "path/to/UMI_count_table.txt"
outdir <- "path/to/output/directory"

# Read count matrix
raw.data <- read.table(dataset, header = TRUE, row.names = 1)
count <- as.matrix(raw.data)

# Run sctransform
data.sctransform <- sctransform::vst(count,n_genes=NULL,min_cells=0)
z = data.sctransform$y

# Compute the mean gene expression from the regularized negative binomial regression
theta = data.sctransform$model_pars_fit[,1]
b0    = data.sctransform$model_pars_fit[,2]
b1    = data.sctransform$model_pars_fit[,3]
median_UMI = median(colSums(count))
mu = exp( b0 + b1*log10(median_UMI) )

# Write output z matrix, the Pearson residuals of the regularized negative binomial regression
out_z.path<- paste(outdir, 'sctransform_z.txt', sep='/')
write.table(z,file=out_z.path,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

# Write output mean gene expression
out_mu.path<- paste(outdir, 'sctransform_mu.txt', sep='/')
write.table(mu,file=out_mu.path,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
