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
vae = VAE(local_csv_dataset.nb_genes, n_batch=local_csv_dataset.n_batches * use_batches)
trainer = UnsupervisedTrainer(vae, local_csv_dataset, train_size=0.75, use_cuda=use_cuda)
trainer.train()
full = trainer.create_posterior(trainer.model, local_csv_dataset, indices=np.arange(len(local_csv_dataset)))
imputed_values = full.sequential().imputation()

# Write output matrix
np.savetxt(outdir + '/scvi_normalization.txt',imputed_values.T,fmt='%.6e',delimiter='\t')
