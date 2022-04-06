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
covariate = args[12]
#adata.X = sparse.csr_matrix(adata.X)
#adata.layers["counts"] = adata.X.copy()
#adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])

# gtemp save
#adata.write_h5ad('scvi.h5ad')


#########################
# Key step
#############################
# trian model on subset (e.g. just human retina)
# project remaining data onto trained model

##############################
# samples = pd.read_csv('/home/mcgaugheyd/git/scEiaD/data/human_ref_samples.txt', header = None)
samples = pd.read_csv(args[11], header = None)
ref_samples = samples.iloc[:,0].to_list()
# selecting by sapmle_accession
if '_samples' in args[11]:
	ref = np.array([s in ref_samples for s in adata.obs.sample_accession])
# selecting by barcodde
elif '_barcode' in args[11]:
	ref = []
	ref_samples = set(ref_samples)	
	for i in adata.obs.index.values:
		if i in ref_samples:
			ref.append(True)
		else:
			ref.append(False)
	ref = np.array(ref)

adata_ref = adata[ref, var_names[0]].copy()
adata_query = adata[~ref, var_names[0]].copy()

adata_ref.X = sparse.csr_matrix(adata_ref.X)
adata_query.X = sparse.csr_matrix(adata_query.X)

scvi.model.SCVI.setup_anndata(adata_ref, batch_key=covariate)#, continuous_covariate_keys = ['percent.mt'])
scvi.model.SCVI.setup_anndata(adata_query, batch_key=covariate) #, continuous_covariate_keys = ['percent.mt'])
arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
	n_latent = n_latent,
)
adata_ref
vae_ref = scvi.model.SCVI(
    adata_ref,
    **arches_params
)
#vae_ref.train(n_epochs = n_epochs, n_epochs_kl_warmup = None)
vae_ref.train(max_epochs = n_epochs, use_gpu=useCuda)
vae_ref

# scANVI
lvae = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    adata=adata_ref,
    labels_key="CellType",
    unlabeled_category="NA",
)
lvae.train(max_epochs=n_epochs, n_samples_per_label=100)

# save the reference model
dir_path = "scVI_HSdroplet_model/" + str(adata.shape[1]) + "HVG_" + str(n_latent) + "ld/"
if len(args) == 14:
	dir_path = args[13]
	lvae.save(dir_path, overwrite=True)

# pull in scVI latent dim
adata_ref.obsm["X_scVI"] = lvae.get_latent_representation()


# remove celltype in query that are not in ref
ct_query = list(adata_query.obs['CellType'].unique())
ct_ref = list(adata_ref.obs['CellType'].unique())
for x in ct_query:
	if x not in ct_ref:
		print('replace ' + x)
		adata_query.obs['CellType'] = adata_query.obs['CellType'].replace([x],'NA')

adata_query
scvi.model.SCANVI.prepare_query_anndata(adata_query, lvae)
vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    lvae,
)

vae_q.train(max_epochs=n_epochs, use_gpu=useCuda, plan_kwargs=dict(weight_decay=0.0))
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()

adata_full = adata_query.concatenate(adata_ref, batch_key = 'bkey')

adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)


obs=pd.DataFrame(adata_full.obs)
obs.to_csv(args[1] + '.' + rand + '.meta.csv')
scvi_latent=pd.DataFrame(adata_full.obsm['X_scVI'])
scvi_latent.to_csv(args[1] + '.' + rand + '.latent.csv')
