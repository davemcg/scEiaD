import sys
import desc
import numpy as np
import pandas as pd
import scanpy.api as sc

args = sys.argv

adata = sc.read_loom(args[1], sparse = True)

desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
desc.log1p(adata)
adata.raw = adata
desc.scale(adata, zero_center = True, max_value = 3)

adata_train = desc.train(adata, 
				   tol=0.1,
                   save_dir="desc_result", 
				   do_tsne=False, 
                   do_umap=False,
				   use_GPU=True,
				   max_iter = 100,
                   save_encoder_weights=True)

obsm_data=pd.DataFrame(adata.obsm["X_Embeded_z0.8"])
obsm_data.to_csv(args[1] + ".csv", sep=",")
