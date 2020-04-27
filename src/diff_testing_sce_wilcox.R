library(Seurat)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)

load(args[1]) # seurat obj
load(args[2]) # cluster
load(args[3]) # cell type prediction
int_sce <-  as.SingleCellExperiment(integrated_obj)

if (all(colnames(int_sce) == (meta$Barcode))) {
	int_sce$cluster <- meta[,2] %>% pull(1)
	int_sce$subcluster <- meta[,3] %>% pull(1)
} else {
	stop('Cluster Barcodes != SCE barcode order')
}

if (args[4] == 'cluster') {
	group = int_sce$cluster
} else if (args[4] == 'subcluster') {
	group = int_sce$subcluster
} else if (args[4] == 'CellType_predict') {
	if (all(colnames(int_sce) == umap$Barcode)){
		group = umap$CellType_predict
	} else { stop('Prediction Barcodes != SCE barcode order')
	}
}

markers_wilcox <- findMarkers(int_sce, 
				group = group, 
				block = int_sce$batch, 
				test="wilcox", 
				BPPARAM=MulticoreParam(as.integer(args[5])))

save(markers_wilcox, file = args[6]) 
