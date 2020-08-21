import sys
import bbknn as bb
import numpy as np
import pandas as pd
import scanpy as sc

args = sys.argv
latent = int(args[2])
adata = sc.read_loom(args[1], sparse = True)
sc.tl.pca(adata, latent)
bb.bbknn(adata, batch_key = 'batch')


obsm_data=pd.DataFrame(adata.varm['PCs'])
obsm_data.to_csv(args[1] + ".csv", sep=",")
