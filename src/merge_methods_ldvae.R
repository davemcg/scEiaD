print('script start')
args <- commandArgs(trailingOnly = TRUE)
save(args, file = 'testing/rldavse.Rdata')
method = args[1]
conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
cat('\n\n')
cat(conda_dir)
cat('\n')
cat(git_dir)
cat('\n\n')
library(glue)
Sys.setenv('RETICULATE_PYTHON' = glue('{conda_dir}/envs/scVI/bin/python') )
library(reticulate)
library(Seurat)
#library(SeuratDisk)
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)

transform = args[2]
covariate = args[3]
latent = args[4] %>% as.numeric()
output_seu_obj_file <- args[6]
#args[5] <- gsub('counts', 'standard', args[5])

load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}

seu <- load_rdata(args[5])
as_SCE_simple <- function(seu){
  t_assays <-  Seurat::Assays(seu)
  sce <- lapply(t_assays,function(x) GetAssayData(seu, slot = 'counts',assay = x ))
  names(sce) <- t_assays
  sce <- SingleCellExperiment(sce)
  rd <- lapply(Reductions(seu), function(x) Embeddings(seu[[x]]) )
  names(rd) <- Reductions(seu)
  reducedDims(sce) <- rd
  metadata(sce) <- seu[[]]
  return(sce)
}
seu_f <- seu[VariableFeatures(seu),]
intg_sce <- as_SCE_simple(seu_f)
batch <- seu_f$batch %>% as.character
intg_adata <- zellkonverter::SCE2AnnData(intg_sce)
use_python(glue('{conda_dir}/envs/scVI/bin/python'))
source_python( glue('{git_dir}/src/ldvae_wrapper.py'))
pfx <- str_split(output_seu_obj_file, '/') %>% sapply(function(x) x[3]) %>% str_remove_all('preFilter.seuratV3.Rdata')
system('mkdir -p ldvae/')
latent_dims <- ldvae_wrapper(adata = intg_adata, batch=batch, prefix=paste0('ldvae/', pfx))
row.names(latent_dims) <- colnames(seu)
colnames(latent_dims) <- paste0("scVI_", 1:ncol(latent_dims))
seu[["scVI"]] <- CreateDimReducObject(embeddings = latent_dims %>% as.matrix(), key = "scVI_", assay = DefaultAssay(seu))
save(seu, file = output_seu_obj_file, compress = FALSE)
