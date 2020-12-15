#%%

import sys
import matplotlib; matplotlib.use('agg')
import torch
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # scvi.data.dataset
from scvi.model import SCVI
from scvi.data import read_loom, setup_anndata  #LDVAE #gone.. 
import scanpy as sc
import os 
import loompy
os.chdir('/data/swamyvs/scEiaD')
#%%
args = sys.argv
## convert loom to anndata
adata = sc.read_loom(args[1])
adata.layers["counts"] = adata.X.copy()
setup_anndata(adata, layer="counts", batch_key="batch")

n_epochs = int(args[2])
lr = float(args[3])
if args[4] == 'True':
	use_cuda = True
else:
	use_cuda = False

n_hidden = int(args[5])
n_latent = int(args[6])
n_layers = int(args[7])


#vae = VAE(loom_dataset.nb_genes, n_batch=loom_dataset.n_batches, dropout_rate = 0.5, n_hidden = n_hidden, n_latent = n_latent, n_layers = n_layers, dispersion='gene-batch')
#vae = VAE(loom_dataset.nb_genes, n_batch=loom_dataset.n_batches, n_latent = n_latent, dropout_rate = 0.1, dispersion='gene-batch')
vae = SCVI(adata,                     
		n_hidden = n_hidden,
		n_latent = n_latent,
		n_layers = n_layers,
		dispersion='gene-batch',
		use_cuda = True)
vae
vae.train(n_epochs = n_epochs, n_epochs_kl_warmup = 1)

# extract
latent = vae.get_latent_representation()


normalized_values = vae.get_normalized_expression()
with open(args[1] + ".csv", 'wb') as f:
	np.savetxt(f, latent, delimiter=",")
if args[8] == 'IMPUTE':
	out = pd.DataFrame(imputed_values)
	out.columns = loom_dataset.gene_names
	out.index = loom_dataset.CellID
	out.to_hdf(args[1] + '.impute.hdf', key = 'index')
if args[8] == 'NORM':
	out = pd.DataFrame(normalized_values) 
	out.columns = loom_dataset.gene_names
	out.index = loom_dataset.CellID
	out.to_hdf(args[1] + '.norm.hdf', key = 'index')



