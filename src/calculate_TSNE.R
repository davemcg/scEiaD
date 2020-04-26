library(tidyverse)
library(Seurat)
library(future)
#library(scran)
plan(strategy = "multicore", workers = 4)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
max_dims = args[2] %>% as.numeric()
perplexity = args[3] %>% as.numeric()

# integrated seurat obj
load(args[4])

print(paste('Method', method))
print(paste('Perplexity', perplexity))

create_TSNE_and_cluster <- function(integrated_obj, 
                                         max_dims = 20, 
                                         perplexity = 200,
                                         reduction = 'pca',
                                         reduction.name = 'ccaTSNE',
                                         reduction.key = 'ccaTSNE_'){
    # TSNE
	print('TSNE 2D Starting')
	integrated_obj <- RunTSNE(integrated_obj,
                              dims = 1:max_dims,
                              check_duplicates = FALSE,
                              dim.embed = 2,
                              perplexity = perplexity,
                              reduction = reduction,
                              learning_rate = 37500,
                              late_exag_coeff = 6,
                              tsne.method = 'FIt-SNE',
                              fast_tsne_path = '/home/mcgaugheyd/git/FIt-SNE/bin/fast_tsne',
                              initialization = 2,
  							  nthreads = 8,
                              reduction.name = reduction.name,
                              reduction.key = reduction.key)  
	integrated_obj
}


if (method == 'CCA'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'ccaTSNE_'
} else if (method == 'scanorama'){
  reduction <- 'scanorama'
  reduction.key <- 'scanoramaTSNE_'
} else if (method == 'harmony'){
  reduction <- 'harmony'
  reduction.key <- 'harmonyTSNE_'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
  reduction.key <- 'mnnTSNE_'
} else if (method == 'none'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'noneTSNE_'
} else if (method == 'combat'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'combatTSNE_'
} else if (method == 'scVI'){
  reduction <- 'scVI'
  reduction.key <- 'scviTSNE_'
} else {
  print(paste0("Why did you pick ", method, "?"))
}
reduction.name <- gsub('_','', reduction.key)

the_reduction = reduction
integrated_obj <- create_TSNE_and_cluster(integrated_obj = integrated_obj, 
                                               max_dims = max_dims,
                                               reduction = the_reduction,
                                               perplexity = perplexity,
                                               reduction.name = reduction.name,
                                               reduction.key = reduction.key)
save(integrated_obj, file = args[5], compress = FALSE )
