# wrapper for python script to load conda
# yes this is crazy
conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RETICULATE_PYTHON = glue('{conda_dir}/envs/scIB/bin/python') )
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

system(paste(glue('{conda_dir}/envs/scIB/bin/python {git_dir}/src/scIB_stats.py '), paste(args, collapse = ' ')))
