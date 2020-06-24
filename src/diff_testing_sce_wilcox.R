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

umap <- umap %>% mutate(CellType = case_when(!CellType %in% c('Astrocytes', 'Fibroblasts', 'Red Blood Cells', 'RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cells', 'Doublet', 'Doublets') ~ CellType))
if (all(colnames(int_sce) == (meta$Barcode))) {
	int_sce$cluster <- meta[,2] %>% pull(1)
	int_sce$subcluster <- meta[,3] %>% pull(1)
    int_sce$CellType <- umap$CellType  
} else {
	stop('Cluster Barcodes != SCE barcode order')
}

if (args[4] == 'cluster') {
	group = int_sce$cluster
} else if (args[4] == 'CellType_predict') {
	if (all(colnames(int_sce) == umap$Barcode)){
		group = umap$CellType_predict
	} else { stop('Prediction Barcodes != SCE barcode order')
	}
} else if (args[4] == 'CellType') {
	if (all(colnames(int_sce) == umap$Barcode)){
		int_sce <- int_sce[,!is.na(int_sce@colData[,'CellType'])]
		group = int_sce$CellType
	} else { stop('Prediction Barcodes != SCE barcode order')
	}
}
if (args[4] == 'subcluster'){
	marker_list <- list()
	for (i in unique(int_sce$cluster)[1:3]){
		sub_int_sce = int_sce[, int_sce$cluster == i]
		print(i)
   		marker_list[[i]] <- findMarkers(sub_int_sce,    
				group = sub_int_sce$subcluster, 
				block = sub_int_sce$batch, 
				pval.type = 'some',
				test="wilcox", 
				min.prop = 0.8,
				BPPARAM=MulticoreParam(as.integer(args[5])))
	}
	markers_wilcox <- do.call(c, unlist(marker_list))
} else { markers_wilcox <- findMarkers(int_sce, 
				group = group, 
				block = int_sce$batch, 
				pval.type = 'some',
				test="wilcox", 
				min.prop = 0.8,
				BPPARAM=MulticoreParam(as.integer(args[5])))
}


markers_summary <- list()
for (i in names(markers_wilcox)){
  print(i)
  temp <- markers_wilcox[[i]][,4:ncol(markers_wilcox[[i]])] %>% as.matrix()
  count <-  apply(temp, 1, function(x) sum(x > 0.7))
  med_auc <- apply(temp, 1, function(x) median(x))
  mean_auc <- apply(temp, 1, function(x) mean(x))
  gene <- row.names(temp)
  cluster <- i
  markers_summary[[i]] <- cbind(gene, count, med_auc, mean_auc, cluster) %>% as_tibble()
}
markers_summary <- markers_summary %>% bind_rows() %>% mutate(count = as.numeric(count), med_auc = as.numeric(med_auc), mean_auc = as.numeric(mean_auc))

save(markers_wilcox, file = args[6])
save(markers_summary, file = args[7])
