# combine datasets by species

# first with 10X/UMI and mouse
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
# metadata
metadata <- read_tsv(args[1])
species <- args[2]

sc_data <- list()
for (i in seq(1,length(args))){
  file = args[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  sc_data[[sample_accession]] <- Seurat::CreateSeuratObject(counts = res_matrix, project = sample_accession)
  sc_data[[sample_accession]] <- NormalizeData(sc_data[[sample_accession]], verbose = FALSE)
  sc_data[[sample_accession]] <- FindVariableFeatures(sc_data[[sample_accession]], selection.method = 'vst', 
                                                      nfeatures = 2000, verbose = FALSE)
  }

samples <- str_extract(args, '(SRS|iPSC_RPE_scRNA_)\\d+')

sra_metadata_extended %>% 
  filter(sample_accession %in% samples)


ref_samples <- sc_data[]
anchors <- FindIntegrationAnchors(sc_data, dims = 1:30)
