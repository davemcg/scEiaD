#!/data/mcgaugheyd/conda/envs/scVI/bin/python

import sys
import loom
import scarches as sca
import scanpy as sc
import pandas as pd

args = sys.argv

adata = sc.read_loom(args[1])
print(adata.X)

n_epochs = int(args[2])
lr = float(args[3])
if args[4] == 'True':
	use_cuda = True
else:
	use_cuda = False

n_hidden = int(args[5])
n_latent = int(args[6])
n_layers = int(args[7])


print(str(n_epochs) + " epochs")
print(str(lr) + " learning rate")
print(str(n_hidden) + " hidden layers")
print(str(n_latent) + " latent dims")
print(str(n_layers) + " layers")

#vae = VAE(loom_dataset.nb_genes, n_batch=loom_dataset.n_batches, n_latent = n_latent, dropout_rate = 0.1, dispersion='gene-batch')
adata = sca.data.normalize_hvg(adata,
		batch_key='batch',
		n_top_genes=adata.X.shape[1])
network = sca.models.scArches(task_name='scEiaD',
                              x_dimension=adata.shape[1],
                              z_dimension=n_latent,
                              architecture=[n_hidden, n_hidden],
                              gene_names=adata.var_names.tolist(),
                              conditions=adata.obs['batch'].unique().tolist(),
                              alpha=lr,
                              loss_fn='nb',
                              model_path="./models/scArches/",
                              )
network.train(adata,
              condition_key='batch',
              n_epochs=n_epochs,              
              batch_size=256,
              save=False,
              retrain=True)

latent_adata = network.get_latent(adata, 'batch')
latent = pd.DataFrame(latent_adata.X)
latent.index = adata.obs.index

latent.to_csv(args[1] + ".csv")


