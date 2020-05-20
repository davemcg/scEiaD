import sys
import tnn
from tnn.tnn import *
import scanpy as sc

args = sys.argv

adata = sc.read_loom(args[1], sparse = True)
dims = float(args[2])
if type == 'counts':
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	adata.raw = adata
# pca
sc.tl.pca(adata)

if 'onlyWE' in args[1]:
	k = 10
else:
	k = 100

# semi-supervised
adata_tnn = adata.copy()
tnn_model = TNN(k = k, 
        distance = 'pn', 
        batch_size = 64, 
        n_epochs_without_progress = 4, 
        epochs = 50, 
        embedding_dims=dims,
        approx = True)
#if 'onlyDROP' in args[1]:
#	tnn_model.fit(X = adata_tnn, 
#		batch_name = "masking_batch", 
#		celltype_name='celltype', 
#		mask_batch = 'missing')
#else:
tnn_model.fit(X = adata_tnn,  batch_name = "batch") 

latent = tnn_model.transform(adata_tnn)

pd.DataFrame(latent).to_csv(args[1] + '.csv',  sep = ',')	
