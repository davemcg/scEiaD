library(reticulate); use_python('/data/mcgaugheyd/conda/envs/pacmap/bin/python3', required = TRUE)
library(tidyverse)
library(Seurat)

py_pacmap <- import('pacmap')

# load(args[1])
# load('pipeline_data/xgboost_predictions/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-6__epochs-15__dist-0.1__neighbors-50__knn-20__umapPredictions.Rdata')

args = commandArgs(trailingOnly=TRUE)
load(args[1])
#  load('seurat_obj/integrated/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-6__epochs-15__preFilter.seuratV3.Rdata')

meta <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode')
X <- meta %>% select(contains('scVI')) %>% as.matrix()

embedding = py_pacmap$PaCMAP(n_dims=as.integer(2), apply_pca = FALSE)
pacmap = embedding$fit_transform(X)

pacmap <- pacmap %>% as_tibble()
pacmap$Barcode <- meta$Barcode
colnames(pacmap)[1:2] <- c('pacmap_1', 'pacmap_2')
save(pacmap, file = args[2])
