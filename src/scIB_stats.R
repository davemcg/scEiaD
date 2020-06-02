# wrapper for python script to load conda
# yes this is crazy
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RETICULATE_PYTHON = '/data/mcgaugheyd/conda/envs/scIB/bin/python')
library(reticulate)


system(paste('/data/mcgaugheyd/conda/envs/scIB/bin/python /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/scIB_stats.py ', paste(args, collapse = ' ')))
