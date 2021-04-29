### python time
import anndata
import sys
import os
import numpy as np
import pandas as pd
import random
import scanpy as sc
from scipy import sparse
import scvi
import torch
# 8 cores
sc.settings.n_jobs = 8
# set seeds
random.seed(234)
scvi.settings.seed = 234
# take in args
args = sys.argv
print(args)
n_epochs = int(args[3])
org = args[5]

# load query data
adata_query = sc.read_h5ad(args[1])
adata_query.layers["counts"] = adata_query.X.copy()
adata_query.layers["counts"] = sparse.csr_matrix(adata_query.layers["counts"])


# Set scVI model path
scVI_model_dir_path = 'scVIprojectionSO_scEiaD_model/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojectionSO__dims-8/' 
# Read in HVG genes used in scVI model
var_names = pd.read_csv(scVI_model_dir_path + '/var_names.csv', header = None)
# cut down query adata object to use just the var_names used in the scVI model training

if org.lower() == 'mouse':
	adata_query.var_names = adata_query.var['gene_name']
	n_missing_genes = sum(~var_names[0].isin(adata_query.var_names))
	dummy_adata = anndata.AnnData(X=sparse.csr_matrix((adata_query.shape[0], n_missing_genes)))
	dummy_adata.obs_names = adata_query.obs_names
	dummy_adata.var_names = var_names[0][~var_names[0].isin(adata_query.var_names)]
	adata_fixed = anndata.concat([adata_query, dummy_adata], axis=1)
	adata_query_HVG = adata_fixed[:, var_names[0]]

adata_query_HVG.obs['batch'] = 'New Data'

scvi.data.setup_anndata(adata_query_HVG, batch_key="batch")
vae_query = scvi.model.SCVI.load_query_data(
    adata_query_HVG, 
    scVI_model_dir_path
)
# project scVI latent dims from scEiaD onto query data
vae_query.train(max_epochs=n_epochs,  plan_kwargs=dict(weight_decay=0.0))
# get the latent dims into the adata
adata_query_HVG.obsm["X_scVI"] = vae_query.get_latent_representation()

# output latent dims to csv
obs=pd.DataFrame(adata_query_HVG.obs)
obsm=pd.DataFrame(adata_query_HVG.obsm["X_scVI"])
features = list(obsm.columns)
obsm.index = obs.index.values
obsm['Barcode'] = obsm.index
obsm['Age'] = 1000
obsm['organism'] = 'x'
# xgboost ML time
from cell_type_predictor import *


CT_predictions = scEiaD_classifier_predict(inputMatrix=obsm, 
                               labelIdCol='ID', 
                               labelNameCol='CellType',  
                               trainedModelFile= os.getcwd() + '/2021_cell_type_ML_all',
                               featureCols=features,  
                               predProbThresh=float(args[4]))
print(CT_predictions['CellType'].value_counts())
CT_predictions.to_csv(args[2])
