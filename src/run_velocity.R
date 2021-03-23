# conda_env <- '/data/swamyvs/anaconda3/bin/python'
# git_dir <- '/data/swamyvs/scEiad_quant/'
# seurat_obj_file <- '/data/OGVFB_BG/scEiaD/2021_03_18/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojectionSO__dims-8__preFilter.scEiaD.seuratV3.Rdata'
# umap_file = '/data/OGVFB_BG/scEiaD/2021_03_18/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojectionSO__dims-8__preFilter.scEiaD__dist-0.1__neighbors-50.umapFilter.predictions.Rdata'
# unspliced_matrix_file <- "/data/OGVFB_BG/scEiaD/2021_02_05/all_species_full_sparse_unspliced_matrix.Rdata"
# spec = "Homo sapiens"
# outfile = paste0('testing/scEiaD_2021_02_10_scvelo_full_sce_nvargenes-ALL', spec ,'_.Rdata')

args <- commandArgs(trailingOnly = T)
# REMEMBER TO LOAD zlib
save(args, file = 'testing/velo.args')
library(glue)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(reticulate)

load_rdata <- function(x){
    load(x)
    env <- ls.str()
    var <- env[!grepl('^x$', env)]
    stopifnot(length(var) == 1)
    return(get(var))
}

conda_env = args[1]
git_dir = args[2]
seurat_obj_file = args[3]
umap_file = args[4]
unspliced_matrix_file = args[5]
outfile = args[7]
recover_dynamics = args[8]=="TRUE"
if (args[6] == 1){
  spec = "Homo sapiens"
} else if (args[6] == 2){
  spec = "Mus musculus"
}
message(spec)

message('Loading Data...\n')
umap = load_rdata(umap_file)
umap_spec <- filter(umap, !is.na(CellType_predict),organism == spec)
labels <- umap_spec %>% select(barcode=Barcode, CellType_predict)
int_seu <- load_rdata(seurat_obj_file)
if(!any('scviUMAP'%in% Reductions(int_seu))){
  umat <- umap[,c('UMAP_1', 'UMAP_2')] %>% as.matrix
  colnames(umat) <- c('UMAP_1', 'UMAP_2')
  rownames(umat) <- umap$Barcode
  common_cells <- intersect(rownames(umat), colnames(int_seu))
  int_seu <- int_seu[,common_cells]
  int_seu[['scviUMAP']] <- CreateDimReducObject(umat[common_cells,], key = 'scviUMAP')
}

unspliced_counts <- load_rdata(unspliced_matrix_file)


# ## extract counts matrix, and subset to common cells 
spliced_counts <- GetAssayData(int_seu,  slot='counts', assay = 'RNA')
common_cells <- intersect(colnames(spliced_counts), colnames(unspliced_counts)) %>% intersect(umap_spec$Barcode)
common_genes <- intersect(rownames(spliced_counts), rownames(unspliced_counts))

spliced_common <- spliced_counts[common_genes, common_cells]
unspliced_common <- unspliced_counts[common_genes, common_cells]
spl_cs = colSums(spliced_common)
uns_cs <- colSums(unspliced_common)
uns_spl_ratio <- uns_cs / (uns_cs + spl_cs) # can use to optionally filter out cells with very low intron quant
# above_min_splice_ratio = uns_spl_ratio > .01 # at least 5% of total expression is intronic quantification

sce <- SingleCellExperiment(list(spliced = spliced_common,
                                 unspliced = unspliced_common
                                 ))
rd <- lapply(int_seu@reductions, function(x) x@cell.embeddings[common_cells,])
names(rd) <- paste0('X_', names(rd))
reducedDims(sce) <- rd

## all metadata( int seu has some, meta data has some, and they overlap)
ccols <-  colnames(umap_spec) [ colnames(umap_spec)%in% colnames(int_seu@meta.data) &  colnames(umap_spec) != 'Barcode' ] # remove duplicate columns 
metadata <- int_seu@meta.data[common_cells,] %>% rownames_to_column('Barcode') %>%  inner_join(umap_spec %>% select(-all_of(ccols)))

for (col in colnames(metadata)){
  sce[[col]] <- metadata[[col]]
  
}
use_python(conda_env)
source(glue('{git_dir}/src/konverter.R'))
source_python(glue('{git_dir}/src/writeH5ad.py') )
message('Converting to anndata...\n')
adata <- SCE2AnnData(sce)
source_python( glue('{git_dir}/src/scvelo_wrapper.py'))

message('Running scvelo...\n')
velo_adata = run_velocity(adata=adata,embedding_key = 'scVI', dkey='scviUMAP', vrank_gb='CellType_predict',n_top_genes = -1, recover_dynamics=recover_dynamics) #  set n_top_genes to -1 for all genes 

message('Converting to SCE...\n')
velo_sce= AnnData2SCE(velo_adata)
save(velo_sce, file = outfile )

