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
#ref = np.array([s in ['E-MTAB-7316_10xv2_Donor1','E-MTAB-7316_10xv2_Donor2','E-MTAB-7316_10xv2_Donor3','OGVFB_Hufnagel_iPSC_RPE_10xv2_None','SRP106476_SMARTerSeq_v3_NA','SRP125998_SMARTSeq_v2_NA','SRP136739_SMARTSeq_v4_NA','SRP151023_10xv2_NA','SRP159286_SCRBSeq_NA','SRP161678_SMARTSeq_v4_NA','SRP170038_SMARTSeq_v2_NA','SRP170761_10xv2_NA','SRP194595_10xv3_Donor1','SRP194595_10xv3_Donor2','SRP194595_10xv3_Donor3','SRP218652_10xv3_donor1','SRP218652_10xv3_donor2','SRP218652_10xv3_donor3','SRP218652_10xv3_donor4','SRP218652_10xv3_donor5','SRP218652_10xv3_donor6','SRP218652_10xv3_donor7','SRP222001_10xv2_retina1','SRP222001_10xv2_retina2','SRP222001_10xv2_retina3','SRP222958_DropSeq_retina2','SRP222958_DropSeq_retina6','SRP222958_DropSeq_retina8','SRP223254_10xv2_NA','SRP223254_10xv2_rep2','SRP238587_10xv2_NA','SRP255195_10xv2_H1','SRP255195_10xv2_H2','SRP255195_10xv2_H3','SRP255195_10xv2_H4','SRP255195_10xv2_H5','SRP255195_10xv3_H1','SRP257883_10xv3_donor_22','SRP257883_10xv3_donor_23','SRP257883_10xv3_donor_24','SRP257883_10xv3_donor_25'] for s in adata.obs.sample_accession])

adata_ref = adata[ref, var_names[0]].copy()
adata_query = adata[~ref, var_names[0]].copy()

adata_ref.X = sparse.csr_matrix(adata_ref.X)
adata_query.X = sparse.csr_matrix(adata_query.X)

scvi.data.setup_anndata(adata_ref, batch_key="batch", continuous_covariate_keys = ['percent.mt'])
scvi.data.setup_anndata(adata_query, batch_key="batch", continuous_covariate_keys = ['percent.mt'])
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


# save the reference model
dir_path = "scVI_HSdroplet_model/" + str(adata.shape[1]) + "HVG_" + str(n_latent) + "ld/"
if len(args) == 13:
	dir_path = args[12]
vae_ref.save(dir_path, overwrite=True)
adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()

adata_query
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    vae_ref,
)

vae_q.train(max_epochs=n_epochs, use_gpu=useCuda, plan_kwargs=dict(weight_decay=0.0))
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()

adata_full = adata_query.concatenate(adata_ref, batch_key = 'bkey')

adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)


obs=pd.DataFrame(adata_full.obs)
obs.to_csv(args[1] + '.' + rand + '.meta.csv')
scvi_latent=pd.DataFrame(adata_full.obsm['X_scVI'])
scvi_latent.to_csv(args[1] + '.' + rand + '.latent.csv')
