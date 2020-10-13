print('script start')
args <- commandArgs(trailingOnly = TRUE)
method = args[1]
# #/**************
#method = 'scVI'
# Sys.setenv('SCIAD_CONDA_DIR'='/data/swamyvs/anaconda3/','SCIAD_GIT_DIR'=getwd() )
# load_rdata <- function(x){
#   load(x)
#   env <- ls.str()
#   var <- env[!grepl('^x$', env)]
#   stopifnot(length(var) == 1)
#   return(get(var))
# }
# seurat_obj <- load_rdata('testing/scvi_test/seurat_obj_no_integration_scvi_small_params_humanonly.Rdata')
# transform = 'counts'
# covariate='batch'
# latent = 10
# #************************/

conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
cat('\n\n')
cat(conda_dir)
cat('\n')
cat(git_dir)
cat('\n\n')


library(glue)
library(loomR)# MUST USE DEVELOPMENT VERSION OF loomR
library(SeuratWrappers)
library(harmony)
library(batchelor)
library(sva)
library(Matrix)
library(data.table)
library(tidyverse)
library(Seurat)
print('library loaded')

transform = args[2]
covariate = args[3]
latent = args[4] %>% as.numeric()
#args[5] <- gsub('counts', 'standard', args[5])
load(args[5])
load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}


run_integration <- function(seurat_obj, method, covariate = 'study_accession', transform = 'standard', latent = 50, file = args[5]){
  # covariate MUST MATCH what was used in build_seurat_obj.R
  # otherwise weird-ness may happen
  # the scaling happens at this level
  # e.g. DO NOT use 'batch' in build_seurat_obj.R then 'study_accession' here
  # scVI ----
  assay <- 'RNA'
  vfeatures <- grep('^MT-', seurat_obj@assays$RNA@var.features, invert =TRUE, value = TRUE)
  if (transform == 'counts'){
    matrix = seurat_obj@assays$RNA@counts[vfeatures, ]
  } else if (transform == 'SCT'){
    assay <- 'SCT'
    vfeatures <- grep('^MT-', seurat_obj@assays$SCT@var.features, invert =TRUE, value = TRUE)
    matrix = seurat_obj@assays$SCT@scale.data[vfeatures, ] %>% Matrix(., sparse = TRUE)
  } else {
    matrix = seurat_obj@assays$RNA@scale.data[vfeatures, ]
  }
  rand <- sample(1e7:9e7,1)
  out <- paste0('scvi_tmp/', method, '_', covariate, '_', transform, '_', length(vfeatures), '_', rand, '_',  latent, '.loom')
  
  # add count to one cell if all are zero
  vfeature_num <- length(vfeatures)
  one0 <- vector(mode = 'numeric', length = vfeature_num)
  one0[2] <- 1
  if (sum(colSums(matrix)==0) > 0){
    matrix[,colSums(matrix) == 0] <- one0
  }
  
  create(filename= out, 
         overwrite = TRUE,
         data = matrix, 
         cell.attrs = list(batch = seurat_obj@meta.data[,covariate],
                           batch_indices = seurat_obj@meta.data[,covariate] %>% 
                             as.factor() %>% 
                             as.numeric()))
  # connect to new loom file, then disconnect...otherwise python call gets borked for 
  # as we are connected into the file on create
  loom <- connect(out, mode = 'r')
  loom$close_all() 
  n_epochs = 5 # use 1e6/# cells of epochs
  lr = 0.001 
  #use_batches = 'True'
  use_cuda = 'True'
  n_hidden = 128 
  n_latent = latent
  n_layers = 2 
  
  scVI_command = paste(glue('{conda_dir}/envs/scVI/bin/./python {git_dir}/src/run_scvi-tools.py'),
                       out,
                       n_epochs,
                       lr,
                       use_cuda,
                       n_hidden,
                       n_latent,
                       n_layers,
                       'FALSE')
  # run scVI     
  print(scVI_command) 
  system(scVI_command)
  # import reduced dim (latent)
  latent_dims <- read.csv(paste0(out, '.csv'), header = FALSE)
  #normalized_values <- fread(paste0(out, '.normalized.csv'), header = FALSE) %>% 
  #  as.matrix()  
  if (latent_dims[1,1] == 'NaN'){
    print('scVI fail, rerunning with fewer hidden dims')
    scVI_command = paste( glue('{conda_dir}/envs/scVI/bin/./python {git_dir}/src/run_scVI.py'), 
                          out, n_epochs, lr, use_cuda, 64, n_latent, n_layers, 'FALSE')
    # run scVI
    print(scVI_command);  system(scVI_command); latent_dims <- read.csv(paste0(out, '.csv'), header = FALSE)
    if (latent_dims[1,1] == 'NaN'){
      print('scVI fail, rerunning with even fewer hidden dims')
      scVI_command = paste(glue('{conda_dir}/envs/scVI/bin/./python {git_dir}/src/run_scVI.py'),
                           out, n_epochs, lr, use_cuda, 50, n_latent, n_layers, 'FALSE')
      # run scVI
      print(scVI_command);  system(scVI_command); latent_dims <- read.csv(paste0(out, '.csv'), header = FALSE)
      if (latent_dims[1,1] == 'NaN'){print("scVI fail again! One more try!"); stop()}
    }
    #   normalized_values <- fread(paste0(out, '.normalized.csv'), header = FALSE) %>% 
    #		as.matrix() 
  }
  
  row.names(latent_dims) <- colnames(seurat_obj)
  colnames(latent_dims) <- paste0("scVI_", 1:ncol(latent_dims))
  
  seurat_obj[["scVI"]] <- CreateDimReducObject(embeddings = latent_dims %>% as.matrix(), key = "scVI_", assay = DefaultAssay(seurat_obj))
  #seurat_obj <- SetAssayData(object = seurat_obj, slot = 'scale.data', new.data = normalized_values)
  #save(normalized_values, file = gsub('.seuratV3.Rdata', '.scVI_scaled.Rdata', args[6]), compress = FALSE)
  # system(paste0('rm ', out, '.normalized.csv'))
  obj <- seurat_obj 
  obj
}
print('running function')
system('mkdir -p scvi_tmp')
if (transform != 'SCT' & method != 'none'){
  integrated_obj <- run_integration(seurat__standard, method, covariate, transform, latent = latent)
} else if (transform == 'SCT' & method != 'none') {
  seurat_list <- seurat__SCT$seurat_list
  if (length(seurat_list) > 1){
    merged <- merge(x = seurat_list[[1]], y = seurat_list[2:length(x = seurat_list)])
  } else {merged <- seurat_list[[1]]}
  merged@assays$SCT@var.features <- seurat__SCT$study_data_features
  DefaultAssay(merged) <- 'SCT'
  merged <- RunPCA(merged, npcs = 100)
  integrated_obj <- run_integration(merged, method, covariate, transform = 'SCT', latent = latent)
} else if (transform != 'SCT' & method == 'none'){
  integrated_obj <- seurat__standard
} else if (transform == 'SCT' & method == 'none'){
  seurat_list <- seurat__SCT$seurat_list
  if (length(seurat_list) > 1){
    merged <- merge(x = seurat_list[[1]], y = seurat_list[2:length(x = seurat_list)])
  } else {merged <- seurat_list[[1]]}
  merged@assays$SCT@var.features <- seurat__SCT$study_data_features
  DefaultAssay(merged) <- 'SCT'
  merged <- RunPCA(merged, npcs = 100)
  integrated_obj <- merged
}
save(integrated_obj, file = args[6], compress = FALSE)


