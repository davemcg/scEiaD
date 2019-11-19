library(tidyverse)
library(Seurat)
library(future)
library(scran)
plan(strategy = "multicore", workers = 4)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
max_dims = args[2] %>% as.numeric()
dist = args[3] %>% as.numeric()
neighbors = args[4] %>% as.numeric()
# integrated seurat obj
load(args[5])

print(paste('Method', method))
print(paste('Max Dims', max_dims))
print(paste('Min Dist', dist))
print(paste('N Neighbors', neighbors))

create_umap_and_cluster <- function(integrated_obj, 
                                    max_dims = 20, 
                                    min_dist = 0.3,
							        n.neighbors = 30,
                                    reduction = 'pca',
                                    reduction.name = 'ccaUMAP',
                                    reduction.key = 'ccaUMAP_',
                                    resolution = 4,
                                    cluster = TRUE){

  # UMAP
  print("UMAP 3D Starting")
  integrated_obj <- RunUMAP(integrated_obj, 
                            dims = 1:max_dims, 
                            min.dist = min_dist,
                            n.components = 3, 
                            n.neighbors = n.neighbors,
                            reduction = reduction, 
                            reduction.name = paste0(reduction.name, '3D'),
                            reduction.key = gsub('_','3D_', reduction.key))
  print("UMAP 2D Starting")
  integrated_obj <- RunUMAP(integrated_obj, 
                            dims = 1:max_dims, 
                            min.dist = min_dist,
                            n.components = 2, 
                            n.neighbors = n.neighbors,
                            reduction = reduction, 
                            reduction.name = reduction.name,
                            reduction.key = reduction.key)
  # clustering 
  # optional cluster on UMAP3D space
  # WAAAAY faster and clusters directly on what you see
  # but also potentially not so valid
  # if you use, drop resolution MUCH lower 0.6-1 or so
  if (!cluster){
    print("Find Neighbors starting")
    #reduction = paste0(reduction.name, '3D')
    #max_dims = 3
    # switch to scran / igraph method as it is about 5X faster
    graph <- buildSNNGraph(Embeddings(integrated_obj, reduction = reduction), transposed = TRUE, d = max_dims)
    cluster <- igraph::cluster_walktrap(graph)$membership
    integrated_obj@meta.data$cluster <- cluster
} else {
	print("Skip clustering")
    integrated_obj@meta.data$cluster <- NA
}
	
  integrated_obj
}


if (method == 'CCA'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
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
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'noneUMAP_'
} else if (method == 'combat'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'combatUMAP_'
} else if (method == 'scVI'){
  reduction <- 'scVI'
  reduction.key <- 'scviUMAP_'
} else {
  print(paste0("Why did you pick ", method, "?"))
}
reduction.name <- gsub('_','', reduction.key)
integrated_obj <- create_umap_and_cluster(integrated_obj = integrated_obj, 
                                          max_dims = max_dims,
                                          reduction = reduction,
                                          min_dist = dist,
										  n.neighbors = neighbors,
                                          reduction.name = reduction.name,
                                          reduction.key = reduction.key,
                                          resolution = 1.5,
                                          cluster = FALSE)

save(integrated_obj, file = args[6], compress = FALSE )

