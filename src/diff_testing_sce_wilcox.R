library(Seurat)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)

load(args[1]) # seurat obj
load(args[2]) # cluster

int_sce <-  as.SingleCellExperiment(integrated_obj)

if (all(colnames(int_sce) == (meta$Barcode))) {
	int_sce$cluster <- meta[,2] %>% pull(1)
	int_sce$subcluster <- meta[,3] %>% pull(1)
} else {
	stop('Cluster Barcodes != SCE barcode order')
}

if (args[3] == 'cluster') {
	group = int_sce$cluster
} else if (args[3] == 'subcluster') {
	group = int_sce$subcluster
}

markers_wilcox <- findMarkers(int_sce, 
				group = group, 
				block = int_sce$batch, 
				test="wilcox", 
				BPPARAM=MulticoreParam(as.integer(args[4])))

save(markers_wilcox, file = args[5]) 
