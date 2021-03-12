args <- commandArgs(trailingOnly = TRUE)
save(args, file = 'testing/cluster.args.Rdata')
library(tidyverse)
library(Seurat)
library(future)
library(scran)
library(glue)
plan(strategy = "multicore", workers = 4)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)

conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
max_dims = args[2] %>% as.numeric()
dist = args[3] %>% as.numeric()
neighbors = args[4] %>% as.numeric()
knn = args[5] %>% as.numeric()
if (args[6] == 'TRUE'){
	cluster = TRUE
} else {cluster = FALSE}
if (args[7] == 'TRUE'){
    umap = TRUE
} else {umap = FALSE}




# integrated seurat obj
load(args[8])

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
                                    knn = 4,
                                    cluster = FALSE,
                                    umap = FALSE){
  set.seed(354)
  if (umap) {
  # UMAP
  #print("UMAP 3D Starting")
  #integrated_obj <- RunUMAP(integrated_obj, 
  #                          dims = 1:max_dims, 
  #                          min.dist = min_dist,
  #                          n.components = 3, 
  #                          n.neighbors = n.neighbors,
  #                          reduction = reduction, 
  #                          reduction.name = paste0(reduction.name, '3D'),
  #                          reduction.key = gsub('_','3D_', reduction.key))
  print("UMAP 2D Starting")
  integrated_obj <- RunUMAP(integrated_obj, 
                            dims = 1:max_dims, 
                            min.dist = min_dist,
                            n.components = 2, 
                            n.neighbors = n.neighbors,
                            reduction = reduction, 
                            reduction.name = reduction.name,
                            reduction.key = reduction.key)
  } else {print("Skip UMAP")}
  # clustering 
  # optional cluster on UMAP3D space
  # WAAAAY faster and clusters directly on what you see
  # but also potentially not so valid
  # if you use, drop knn MUCH lower 0.6-1 or so
  if (cluster){
	if (knn >= 2){
	
   		 print("Clustering begining!")
   	 	#reduction = paste0(reduction.name, '3D')
   	 	#max_dims = 3
   	 	# switch to scran / igraph method as it is about 5X faster
		for (k in knn){
  	 	  graph <- buildSNNGraph(Embeddings(integrated_obj, reduction = reduction), 
		  transposed = TRUE, k = k, d = max_dims, type = 'jaccard')
   	  	  cluster <- igraph::cluster_louvain(graph)$membership
   	 	  integrated_obj@meta.data[,paste0('cluster_knn', k)] <- cluster
		  integrated_obj@misc[[paste0("Graph_knn", k)]] <- graph

		  # do subclustering on each cluster
		  
	  integrated_obj@meta.data[,paste0('subcluster_knn', k)] <- "MSG"
	}
	} else {
		# run parc
		out_embeddings_file = paste(args[8], method, max_dims, dist, neighbors, knn, '.csv', sep = '_') 
		in_embeddings_file = paste0(out_embeddings_file, 'RUN')
		write.csv(Embeddings(integrated_obj, reduction = reduction), file = out_embeddings_file)
		system(paste( glue('{conda_dir}/envs/parc/bin/python {git_dir}/src/run_parc.py'), out_embeddings_file, knn, in_embeddings_file))
		clusters = read.csv(in_embeddings_file)
		met <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
        met_clu <- left_join(met, clusters, by = 'Barcode')
		integrated_obj@meta.data[,paste0('cluster_knn', knn)] <- met_clu$X0
		integrated_obj@misc[[paste0("Graph_knn", knn)]] <- NA
		parc_k = knn + 2
		integrated_obj@meta.data[,paste0('subcluster_knn', parc_k)] <- "MSG"
		system(paste0('rm ', out_embeddings_file))
		system(paste0('rm ', in_embeddings_file))
	}
} else {
    #integrated_obj@meta.data$cluster <- NA
}
  integrated_obj
}


if (method == 'CCA'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'ccaUMAP_'
} else if (method == 'liger'){
  reduction <- 'iNMF'
  reduction.key <- 'iNMFUMAP_'
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
  if (length(integrated_obj@reductions) == 0) {
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'noneUMAP_'
} else if (method == 'combat'){
  if (length(integrated_obj@reductions) == 0) {
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'combatUMAP_'
} else if (method == 'liger'){
  reduction <- 'iNMF'
  reduction.key <- 'ligerUMAP_'
} else if (grepl('scVI', method)){
  reduction <- 'scVI'
  reduction.key <- 'scviUMAP_'
} else if (method == 'magic') {
  reduction <- 'magic'
  reduction.key <- 'magicUMAP_'
} else if (method == 'desc') {
  reduction <- 'desc'
  reduction.key <- 'descUMAP_'
} else {
  reduction <- method
  reduction.key <- paste0(method,'UMAP_')
}
reduction.name <- gsub('_','', reduction.key)
integrated_obj <- create_umap_and_cluster(integrated_obj = integrated_obj, 
                                          max_dims = max_dims,
                                          reduction = reduction,
                                          min_dist = dist,
										  n.neighbors = neighbors,
                                          reduction.name = reduction.name,
                                          reduction.key = reduction.key,
                                          knn = knn,
                                          cluster = cluster,
										  umap = umap)

integrated_obj@assays$RNA@data <- as.matrix(c(0))
integrated_obj@assays$RNA@scale.data <- as.matrix(c(0))
save(integrated_obj, file = args[9], compress = TRUE)

