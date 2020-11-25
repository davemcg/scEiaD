import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import scvi
import skmisc
    
def ldvae_wrapper(adata, batch, prefix):
    adata.obs['batch']=batch
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10e4)
    sc.pp.log1p(adata)
    #adata.raw = adata # freeze the state in `.raw`
    scvi.data.setup_anndata(adata, layer="counts", batch_key = 'batch')
    model = scvi.model.LinearSCVI(adata, n_latent=10)
    model.train(n_epochs=250, lr = 5e-3, frequency = 10)
    training_info = pd.DataFrame.from_dict(
                        {'train_elbo' : model.trainer.history['elbo_train_set'][1:],
                        'test_elbo' : model.trainer.history['elbo_test_set'][1:]})
    embedding = pd.DataFrame(model.get_latent_representation()).rename(columns= lambda x: 'Z_'+str(x))
    loadings = model.get_loadings()
    embedding.to_csv(f'{prefix}_embeddings.csv.gz', index=False)
    loadings.to_csv(f'{prefix}_loadings.csv.gz')
    training_info.to_csv(f'{prefix}_traininginfo.csv.gz',  index=False)
    return embedding