library(tidyverse)
library(Seurat)
library(future)
plan(strategy = "multicore", workers = 6)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
max_dims = args[2] %>% as.numeric()
# cell labels
cell_info_labels <- read_tsv(args[3])
# integrated_obj
load(args[4])


create_umap_and_cluster <- function(integrated_obj, 
                                    max_dims = 20, 
                                    reduction = 'pca',
                                    reduction.name = 'ccaUMAP',
                                    reduction.key = 'ccaUMAP_'){
  # UMAP
  integrated_obj <- RunUMAP(integrated_obj, 
                            dims = 1:max_dims, 
                            reduction = reduction, 
                            reduction.name = reduction.name,
                            reduction.key = reduction.key)
  # clustering 
  integrated_obj <- FindNeighbors(integrated_obj, 
                                  reduction = reduction,
                                  dims = 1:max_dims, 
                                  nn.eps = 0.5)
  integrated_obj <- FindClusters(integrated_obj, 
                                 #resolution = c(0.6,0.8,1,3),
                                 resolution = 1,
                                 save.SNN = TRUE,
                                 do.sparse = TRUE,
                                 algorithm = 2,
                                 random.seed = 23)
  integrated_obj
}

if (method == 'CCA'){
  reduction <- 'pca'
  reduction.key <- 'ccaUMAP_'
} else if (method == 'scanorama'){
  reduction <- 'scanorama'
  reduction.key <- 'scanoramaUMAP_'
} else if (method == 'harmony'){
  reduction <- 'harmony'
  reduction.key <- 'harmonyUMAP_'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
  reduction.key <- 'mnnUMAP_'
} else if (method == 'none'){
  reduction <- 'pca'
  reduction.key <- 'noneUMAP_'
} else if (method == 'combat'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'combatUMAP_'
} else {
  print(paste0("Why did you pick ", method, "?"))
}
reduction.name <- gsub('_','', reduction.key)
integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          max_dims,
                                          reduction,
                                          reduction.name = reduction.name,
                                          reduction.key = reduction.key)

save(integrated_obj, file = args[5], compress = FALSE )

