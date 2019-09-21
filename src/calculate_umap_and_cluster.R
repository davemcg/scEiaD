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
                                 #resolution = c(0.1,0.3,0.6,0.8,1,2,3,4,5),
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
} else {
  print(paste0("Why did you pick ", method, "?"))
}
reduction.name <- gsub('_','', reduction.key)
integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          20,
                                          reduction,
                                          reduction.name = reduction.name,
                                          reduction.key = reduction.key)

save(integrated_obj, file = args[4], compress = FALSE )

# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(integrated_obj[[reduction.name]]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., cell_info_labels %>% dplyr::rename(Barcode = value))
umap$Method <- method
colnames(umap)[2:3] <- c('UMAP_1', 'UMAP_2')

save(umap, file = args[5])
