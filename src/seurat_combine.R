# combine datasets by species

# first with 10X/UMI and mouse
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
options(future.globals.maxSize = 40000 * 1024^2)
# metadata
metadata <- read_tsv(args[1])
species <- args[2]

sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  sc_data[[sample_accession]] <- seurat_obj
}

sc_data_features <- SelectIntegrationFeatures(object.list = sc_data, nfeatures = 3000, verbose = FALSE)

sc_data <- PrepSCTIntegration(object.list = sc_data, anchor.features = sc_data_features, verbose = FALSE)

sc_data_anchors <- FindIntegrationAnchors(object.list = sc_data, normalization.method = 'SCT', anchor.features = sc_data_features, verbose = FALSE)

sc_data_integrated <- IntegrateData(anchorset = sc_data_anchors, normalization.method = 'SCT', verbose = FALSE)
