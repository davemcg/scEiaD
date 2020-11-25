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
meta = args[2]
ld = args[3]
scVI_model_dir_path = args[4]

var_names = pd.read_csv(scVI_model_dir_path + '/var_names.csv', header = None)
var_names[0]

# load anndata
adata = sc.read_h5ad(args[1])

# create new anndata object with the var_names in the same order as the model
adata_subGenes = adata[:, var_names[0]]

vae_q = scvi.model.SCVI.load_query_data(
    adata_subGenes,
    scVI_model_dir_path,
)

# forcing adata.X to be CSR
adata_subGenes.X = sp.sparse.csr_matrix(adata_subGenes.X)
vae_q.train(n_epochs=5, weight_decay=0.0, n_epochs_kl_warmup=1, lr = 0.001)

# extract LD from model
adata_subGenes.obsm["X_scVI"] = vae_q.get_latent_representation()

# write meta and LD
# to jam back into the seurat obj
obs=pd.DataFrame(adata_subGenes.obs)
obs.to_csv(meta)
scvi_latent=pd.DataFrame(adata_subGenes.obsm['X_scVI'])
scvi_latent.to_csv(ld)
