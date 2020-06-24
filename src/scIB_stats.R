# wrapper for python script to load conda
# yes this is crazy
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RETICULATE_PYTHON = '/data/mcgaugheyd/conda/envs/scIB/bin/python')
library(reticulate)

method = args[1]
if (method == 'CCA'){
  reduction <- 'pca'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
} else if (method == 'none'){
  reduction <- 'pca'
} else if (method == 'combat'){
  reduction <- 'pca'
} else if (method == 'liger'){
  reduction <- 'iNMF'
} else if (method == 'scVI'){
   reduction <- 'scVI'
} else {
  print("GUESSING!")
  reduction <- method
}

args[1] <- paste0('X_', tolower(reduction))

system(paste('/data/mcgaugheyd/conda/envs/scIB/bin/python /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/scIB_stats.py ', paste(args, collapse = ' ')))
