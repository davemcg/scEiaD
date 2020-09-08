# Use >= R/4.0
library(tidyverse)
library(scran)
library(scuttle)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
load(args[1]) 
load(args[2]) # meta, cluster info
load(args[3]) # umap filter predictions
method = args[4]
out <- args[5]

integrated_obj@meta.data$cluster <- meta[,2] %>% pull() %>% as.factor()
colnames(meta)[c(2,3)] <- c('cluster','subcluster')
umap <- umap %>% select(-cluster, -subcluster) %>% left_join(meta, by = 'Barcode')
seurat <- integrated_obj[, umap$Barcode]




quick_label_cluster <- function(umap){
	quick_label <-	umap %>%
	  #mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>%
	  #mutate(CellType = CellType_predict) %>%
 	 group_by(cluster, CellType) %>%
 	 filter(!is.na(CellType), !is.na(cluster)) %>%
 	 summarise(Count = n(), x = mean(UMAP_1), y = mean(UMAP_2),
            Organism = list(unique(organism)),
            study_accession = list(unique(study_accession))) %>%
	  mutate(freq = Count / sum(Count)) %>%
 	 filter(freq > 0.25) %>%
  	ungroup() %>%
 	 group_by(cluster) %>%
  	top_n(3, -freq) %>%
 	 #filter(Count > 100) %>%
  	summarise(seurat_cluster_CellType = paste0(cluster, ': ', paste0(CellType, collapse = ' ')),
            x = mean(x), y = mean(y),
            Count = sum(Count)) %>%
	filter(grepl('Precu|RPC|Ganglia|Amacrine|Bipolar|Rods|Cones|Muller|Neurog|Horizon', seurat_cluster_CellType)) %>%
	unique()
	
}


# cut down umap and seurat obj to only labelled cells in major retina classes
umapF <- umap %>% filter(!is.na(CellType)) %>% filter(grepl('Precu|RPC|Ganglia|Amacrine|Bipolar|Rods|Cones|Muller|Neurog|Horizon', CellType))
sCT = seurat[, umapF$Barcode]
sCT[["UMAP"]] <- CreateDimReducObject(embeddings =
								umapF[,2:3] %>% as.matrix(),
                                key = "UMAP_",
                                assay = DefaultAssay(sCT))

# build labelled cluster obj
cl_quick <- quick_label_cluster(umapF)
umapQL <- umap %>% filter(cluster %in% cl_quick$cluster) %>% left_join(cl_quick, by = 'cluster')
sCT_CL = seurat[, umapQL$Barcode]
sCT_CL[["UMAP"]] <- CreateDimReducObject(embeddings =
                                umapQL[,2:3] %>% as.matrix(),
                                key = "UMAP_",
                                assay = DefaultAssay(sCT_CL))


run_tscan <- function(seurat, group, reduction = 'scVI'){
	sce <- as.SingleCellExperiment(seurat)
	sce$group <- group

	agg <- aggregateAcrossCells(sce, ids=sce$group)
	centroids <- reducedDim(agg, toupper(reduction))
	centroids_UMAP <- reducedDim(agg, "UMAP")

	mst <- createClusterMST(centroids,  outgroup = FALSE )	
	edges <- connectClusterMST(centroids_UMAP, mst, combined=FALSE)
	start_options <- names(igraph::V(mst)[igraph::degree(mst)==1])
	start = (grep('RPC|Neuro', start_options, value = TRUE) %>% factor(., levels = c('Early RPCs','RPCs','Late RPCs', 'Neurogenic')))[1]
	if (is.na(start)) {
		start <- (grep('RPC|Neuro', start_options, value = TRUE) )[1]
	}
	if (is.na(start) || length(start) == 0){	
		start <- names(igraph::V(mst)[igraph::degree(mst)==1])[1]
	} 
	ordering <- orderClusterMST(reducedDim(sce, type = toupper(reduction)), sce$group, centroids, mst)
	
	row.names(ordering) <- colnames(sce)
	out <- list(mst = mst, edges = edges, ordering = ordering)
	out
}	

if (method == 'CCA'){
  reduction <- 'pca'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
} else if (method == 'none'){
  reduction <- 'pca'
} else if (method == 'combat'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
} else if (method == 'liger'){
  reduction <- 'iNMF'
} else if (method == 'scVI'){
   reduction <- 'scVI'
} else {
  print("GUESSING!")
  reduction <- method
}

rm(integrated_obj)

mst_cluster <- run_tscan(sCT_CL, umapQL$seurat_cluster_CellType, reduction)
mst_celltype <- run_tscan(sCT, umapF$CellType, reduction)
save(mst_cluster, mst_celltype, umapQL, umapF, file = out)
