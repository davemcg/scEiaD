# R

args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)

# tx file to change gene ID to gene name
tx <- read_tsv(args[1], col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')
# load bustools count data for UMI data
load(args[2])
# sample name
sample_accession <- args[3]
# output Rdata
output_file <- args[4]

# replace ensgene ID with gene names
# macaca already fine
if (!grepl('macaca', args[1], ignore.case = TRUE)){
	row.names(res_matrix) <- row.names(res_matrix) %>% 
		enframe(value = 'id') %>% 
		left_join(., tx, by = 'id') %>% 
		pull(gene)
}

# add sample ID to UMI
colnames(res_matrix) <- paste0(colnames(res_matrix), '_', sample_accession)

seurat_obj <- CreateSeuratObject(counts = res_matrix, project = sample_accession)
seurat_obj <- SCTransform(object = seurat_obj)
save(seurat_obj, file = output_file)
