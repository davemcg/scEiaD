# combine datasets by species

# first with 10X/UMI and mouse
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
plan(strategy = "multicore", workers = 16)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 40000 * 1024^2)

# load well (not UMI/droplet) data
load(args[2])

rdata_files = args[3:length(args)]
sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  sc_data[[sample_accession]] <- seurat_obj
}

sc_data_features <- SelectIntegrationFeatures(object.list = sc_data, nfeatures = 3000, verbose = FALSE)

sc_data <- PrepSCTIntegration(object.list = sc_data, anchor.features = sc_data_features, verbose = FALSE)
# this one takes HOURS
sc_data_anchors <- FindIntegrationAnchors(object.list = sc_data, normalization.method = 'SCT', anchor.features = sc_data_features, verbose = TRUE)

sc_data_integrated <- IntegrateData(anchorset = sc_data_anchors, normalization.method = 'SCT', verbose = FALSE)

# PCA, UMAP
sc_data_integrated <- RunPCA(sc_data_integrated, npcs = 100, ndims.print = 1:5, nfeatures.print = 5))
sc_data_integrated <- RunUMAP(sc_data_integrated, dims = 1:75, min.dist = 0.75)

save(sc_data_integrated, file = args[1])


