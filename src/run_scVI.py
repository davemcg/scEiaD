#!/data/mcgaugheyd/conda/envs/scVI/bin/python

import sys
import matplotlib; matplotlib.use('agg')
import torch
import os
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE, LDVAE
from scvi.inference import UnsupervisedTrainer
from scvi.inference.posterior import Posterior
from scvi.dataset import LoomDataset


args = sys.argv

loom_dataset = LoomDataset(args[1], save_path = os.getcwd())

n_epochs = int(args[2])
lr = float(args[3])
if args[4] == 'True':
	use_cuda = True
else:
	use_cuda = False

n_hidden = int(args[5])
n_latent = int(args[6])
n_layers = int(args[7])

vae = VAE(loom_dataset.nb_genes, n_batch=loom_dataset.n_batches,
        n_hidden = n_hidden, n_latent = n_latent, n_layers = n_layers, dispersion='gene-batch')

trainer = UnsupervisedTrainer(vae, loom_dataset, train_size=1.0, use_cuda=use_cuda)
#n_epochs = 5
trainer.train(n_epochs=n_epochs)
 
# extract
full = trainer.create_posterior(trainer.model, loom_dataset, indices=np.arange(len(loom_dataset)))
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()
 
np.savetxt(args[1] + ".csv", latent, delimiter=",")
 

