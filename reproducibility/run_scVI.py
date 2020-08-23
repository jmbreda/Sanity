

# Import libraries
from scvi.dataset import CsvDataset
import os
import numpy as np
import pandas as pd
from scvi.models import *
from scvi.inference import UnsupervisedTrainer
import torch

# Define input file and output directory
dataset = "path/to/UMI_count_table.csv.gz"
dataset_dir = "path/to/"
outdir = "path/to/output/directory/"

# Read count matrix with all genes (from https://github.com/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb)
local_csv_dataset = CsvDataset(dataset,save_path=dataset_dir,compression='gzip',new_n_genes=False)

# Process data (from https://github.com/YosefLab/scVI/blob/master/tests/notebooks/basic_tutorial.ipynb)
use_batches=False
use_cuda=True
n_epochs=400
vae = VAE(local_csv_dataset.nb_genes, n_batch=local_csv_dataset.n_batches * use_batches)
trainer = UnsupervisedTrainer(vae, local_csv_dataset, train_size=0.75, use_cuda=use_cuda)
trainer.train(n_epochs=n_epochs)
full = trainer.create_posterior(trainer.model, local_csv_dataset, indices=np.arange(len(local_csv_dataset)))
#imputed_values = full.sequential().imputation()
normalized_values = full.sequential().get_sample_scale()

# Write output matrix
np.savetxt(outdir + '/scvi_normalization.txt',normalized_values.T,fmt='%.6e',delimiter='\t')
