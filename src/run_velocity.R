#args <- commandArgs(trailingOnly = T)
setwd('/data/swamyvs/scEiaD/')
conda_dir <- '/data/swamyvs/anaconda3/'
git_dir <- '/data/swamyvs/scEiaD/'
library(glue)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)
seurat_obj_file <- "/data/OGVFB_BG/scEiaD_mcg/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features5000__counts__universe__batch__scVIprojectionSO__dims10__preFilter.scEiaDprojected__mindist0.1__nneighbors500.umap.Rdata"

unspliced_matrix_file <- 'pipeline_data/clean_quant/all_species_full_sparse_unspliced_matrix.Rdata'

load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}

int_seu <- load_rdata(seurat_obj_file)
unspliced_counts <- load_rdata(unspliced_matrix_file)
cluster_assignment <- load_rdata('/data/OGVFB_BG/scEiaD_mcg/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features5000__counts__universe__batch__scVIprojectionSO__dims10__preFilter.scEiaDprojected__knn0.6.cluster.Rdata')
tcol = colnames(cluster_assignment)[grepl("^cluster_knn", colnames(cluster_assignment))]
cluster_assignment= cluster_assignment %>% 
  dplyr::rename(barcode = Barcode, cluster= !!tcol) %>%
  mutate(cluster = paste0('clstr_', cluster))
## extract counts matrix, and 
spliced_counts <- GetAssayData(int_seu,  slot='counts', assay = 'RNA')
common_cells <- intersect(colnames(spliced_counts), colnames(unspliced_counts))
common_genes <- intersect(rownames(spliced_counts), rownames(unspliced_counts))

spliced_common <- spliced_counts[common_genes, common_cells]
unspliced_common <- unspliced_counts[common_genes, common_cells]
spl_cs = colSums(spliced_common)
uns_cs <- colSums(unspliced_common)
uns_spl_ratio <- uns_cs / (uns_cs + spl_cs)
#above_min_splice_ratio = uns_spl_ratio > .01 # at least 5% of total expression is intronic quantification

sce <- SingleCellExperiment(list(spliced = spliced_common,
                                 unspliced = unspliced_common
                                 ))
rd <- lapply(int_seu@reductions, function(x) x@cell.embeddings[common_cells,])
names(rd) <- paste0('X_', names(rd))
reducedDims(sce) <- rd
metadata <- int_seu@meta.data[common_cells,] %>% rownames_to_column('barcode') %>%  inner_join(cluster_assignment)
for (col in colnames(metadata)){
  sce[[col]] <- metadata[[col]]
  
}
adata <- SCE2AnnData(sce)
use_python(glue('{conda_dir}/bin/python'))
source_python( glue('{git_dir}/src/scvelo_wrapper.py'))

velo_adata = run_velocity(adata=adata,embedding_key = 'scVI', dkey='scviUMAP', vrank_gb='cluster')
velo_sce= AnnData2SCE(velo_adata)
save(velo_sce, file = 'testing/TabulaDroplet_vs_scvelo_full_sce.Rdata')



