import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import os

dataset = 'Simulated_Branched_Random_Walk'

My_norm = ['True','RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI','SanityErrorbar'];

# Get initial conf
# Load True distance
D = np.loadtxt('data/Simulated_Branched_Random_Walk_True_Euclidean_dist.txt')

# compute tsne on True distances
X0 = TSNE(n_components=2,perplexity=1200,metric='precomputed').fit_transform(D)

# Compute tsne on distance matrices from different normalization
for my_norm in My_norm:
	file_in = 'data/Simulated_Branched_Random_Walk_' + my_norm + '_Euclidean_dist.txt'
	file_out = 'data/Simulated_Branched_Random_Walk_' + my_norm + '_tsne.txt'

	# load distance matrix
	D = np.loadtxt(file_in)

	# compute tsne
	X = TSNE(n_components=2,perplexity=600,metric='precomputed',init=X_0).fit_transform(D)
	np.savetxt(file_out,X,fmt='%.6f')
