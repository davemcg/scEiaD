library(tidyverse)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
# cell labels
load(args[2])
# integrated_obj
load(args[3])


create_umap_and_cluster <- function(integrated_obj, 
                                  max_dims = 20, 
                                  reduction_name = 'pca', 
                                  reduction_name_key = 'ccaUMAP_'){
  # UMAP
  integrated_obj <- RunUMAP(integrated_obj, 
                            dims = 1:max_dims, 
                            reduction = reduction_name, 
                            reduction.key = reduction_name_key)
  # clustering 
  integrated_obj <- FindNeighbors(integrated_obj, 
                                  dims = 1:max_dims, 
                                  nn.eps = 0.5, 
                                  reduction = reduction_name)
  integrated_obj <- FindClusters(integrated_obj, 
                                 #resolution = c(0.1,0.3,0.6,0.8,1,2,3,4,5),
                                 save.SNN = TRUE,
                                 do.sparse = TRUE,
                                 algorithm = 2,
                                 random.seed = 23)
  integrated_obj
}

if (method == 'CCA'){
  reduction_name_key <- 'ccaUMAP_'
  integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          20,
                                          'pca',
                                          reduction_name_key = reduction_name_key)
} else if (method == 'scanorama'){
  reduction_name_key <- 'scanoramaUMAP_'
  integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          20,
                                          'scanorama',
                                          reduction_name_key = reduction_name_key)
} else if (method == 'harmony'){
  reduction_name_key <- 'harmonyUMAP_'
  integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          20,
                                          'harmony',
                                          reduction_name_key = reduction_name_key)
} else if (method = 'fastMNN'){
  reduction_name_key <- 'mnnUMAP_'
  integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          20,
                                          'mnn',
                                          reduction_name_key = reduction_name_key)
} else {
  print(paste0("Why did you pick ", method, "?"))
}

save(integrated_obj, file = args[4], compress = FALSE )

# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(integrated_obj[[reduction_name_key]]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., meta %>% dplyr::rename(Barcode = value))
umap$Method <- method
colnames(umap)[2:3] <- c('UMAP_1', 'UMAP_2')

save(umap, file = args[5])
