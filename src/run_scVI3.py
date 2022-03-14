import sys
import os
import numpy as np
import pandas as pd
import random
import scanpy as sc
from scipy import sparse
import scvi
import torch

sc.settings.n_jobs = 8
random.seed(234)
scvi.settings.seed = 234

args = sys.argv
print(len(args))
print(args)
n_epochs = int(args[4])
lr = float(args[5])
if args[6] == 'True':
	useCuda = True
else:
	useCuda = False

n_hidden = int(args[7])
n_latent = int(args[8])
n_layers = int(args[9])

var_names = pd.read_csv(args[2], header = None)
rand = args[3]
adata = sc.read_h5ad(args[1])
print('adata loaded')
##############################
# samples = pd.read_csv('/home/mcgaugheyd/git/scEiaD/data/human_ref_samples.txt', header = None)

adata_ref = adata[:, var_names[0]].copy()
print('adata ref made')
adata_ref.X = sparse.csr_matrix(adata_ref.X)
print('adata converted to csr')
scvi.data.setup_anndata(adata_ref, batch_key="batch", continuous_covariate_keys = ['percent.mt'])
scvi_params = dict(
    n_layers=2,
    dispersion = 'gene',
    n_latent = n_latent,
)


adata_ref
vae_ref = scvi.model.SCVI(
    adata_ref,
    **scvi_params
)
#vae_ref.train(n_epochs = n_epochs, n_epochs_kl_warmup = None)
vae_ref.train(max_epochs = n_epochs, use_gpu=useCuda)
vae_ref


# save the reference model
dir_path = args[11]
vae_ref.save(dir_path, overwrite=True)
adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()

# 
obs=pd.DataFrame(adata_ref.obs)
obs.to_csv(args[1] + '.' + rand + '.meta.csv')
scvi_latent=pd.DataFrame(adata_ref.obsm['X_scVI'])
scvi_latent.to_csv(args[1] + '.' + rand + '.latent.csv')
