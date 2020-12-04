import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scipy as sp
import scanpy as sc
import scvi
import random
import sys
random.seed(234)
sc.settings.n_jobs = 8
args = sys.argv
meta = args[3]
ld = args[4]
scVI_model_dir_path = args[5]

var_names = pd.read_csv(scVI_model_dir_path + '/var_names.csv', header = None)
var_names[0]

# load anndata
adata_ref = sc.read_h5ad(args[1])
adata_query = sc.read_h5ad(args[2])
# create new anndata object with the var_names in the same order as the model
adata_subGenes_ref = adata_ref[:, var_names[0]]
adata_subGenes_query = adata_query[:, var_names[0]]
vae_q = scvi.model.SCVI.load_query_data(
    adata_subGenes_query,
    scVI_model_dir_path,
)

# forcing adata.X to be CSR
adata_subGenes_query.X = sp.sparse.csr_matrix(adata_subGenes_query.X)
adata_subGenes_ref.X = sp.sparse.csr_matrix(adata_subGenes_ref.X)
# train
vae_q.train(n_epochs=5, weight_decay=0.0, n_epochs_kl_warmup=1, lr = 0.001)

# concentate with existing
adata_full = adata_subGenes_query.concatenate(adata_subGenes_ref)
adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)
# extract LD from model
adata_full.obsm["X_scVI"] = vae_q.get_latent_representation()

# write meta and LD
# to jam back into the seurat obj
obs=pd.DataFrame(adata_full.obs)
obs.to_csv(meta)
scvi_latent=pd.DataFrame(adata_full.obsm['X_scVI'])
scvi_latent.to_csv(ld)
