library(Seurat)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(BiocParallel)
git_dir = Sys.getenv('SCIAD_GIT_DIR')
args = commandArgs(trailingOnly=TRUE)

load(args[1]) # seurat obj
load(args[2]) # cluster
load(args[3]) # cell type prediction
chunk_to_run <- as.integer(args[8])
how_many_chunks <- as.integer(args[9])
integrated_obj <- NormalizeData(integrated_obj)
x <- CreateSeuratObject(integrated_obj@assays$RNA@data, assay = 'RNA')
integrated_obj@assays$RNA <- x@assays$RNA
rm(x)
# remove non-tissue from diff testing
umap <- umap %>% filter(Source == 'Tissue')
# hand fix some labels
source(glue::glue('{git_dir}/src/tweak_celltype_labels.R'))
umap <- hand_fixer(umap)


int_sce <-  as.SingleCellExperiment(integrated_obj)
int_sce <- int_sce[,umap$Barcode]

meta <- umap %>% select(Barcode) %>% left_join(meta)


# save mem
rm(integrated_obj)

# make chunk
chunk_coords <- list()
for (i in seq(1:how_many_chunks)){
    chunk <- round(nrow(int_sce)/how_many_chunks)
    end = chunk * i
    start = (chunk * i) - chunk + 1
    if (i == how_many_chunks){
        end = nrow(int_sce)
    }
    chunk_coords[[i]] <- c(start,end)
}
chunk_coords <- do.call(rbind, chunk_coords) %>% as.matrix()
int_sce <-int_sce[seq(chunk_coords[chunk_to_run,][1], chunk_coords[chunk_to_run,][2]),]
print("GO DIFF")
# run diff testing
if (all(colnames(int_sce) == (meta$Barcode))) {
	int_sce$cluster <- meta[,2] %>% pull(1)
	#int_sce$subcluster <- meta[,3] %>% pull(1)
    int_sce$CellType_predict <- umap$CellType_predict
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
				block = sub_int_sce$study_accession, 
				pval.type = 'some',
				test="wilcox", 
				min.prop = 0.5,
				BPPARAM=MulticoreParam(as.integer(args[5])))
	}
	markers_wilcox <- do.call(c, unlist(marker_list))
} else { markers_wilcox <- findMarkers(int_sce, 
				group = group, 
				block = int_sce$study_accession, 
				pval.type = 'some',
				test="wilcox",
				min.prop = 0.3,
				BPPARAM=MulticoreParam(as.integer(args[5])))
		markers_t <- findMarkers(int_sce,
                group = group,
                block = int_sce$study_accession,
                pval.type = 'some',
                test="t",
                min.prop = 0.3,
                BPPARAM=MulticoreParam(as.integer(args[5])))
}



markers_summary <- list()
for (i in names(markers_wilcox)){
    print(i)
    temp <- markers_wilcox[[i]][,4:ncol(markers_wilcox[[i]])] %>% as.matrix()
    count <-  apply(temp, 1, function(x) sum(x > 0.7, na.rm = TRUE))
	pval = markers_wilcox[[i]][,1]
	FDR = markers_wilcox[[i]][,2]
    med_auc <- apply(temp, 1, function(x) median(x, na.rm = TRUE))
    mean_auc <- apply(temp, 1, function(x) mean(x, na.rm = TRUE))
    gene <- row.names(temp)
    cluster <- i
    markers_summary[[i]] <- cbind(gene, pval, FDR, count, med_auc, mean_auc, cluster) %>% as_tibble()
}

markers_summary_t <- list()
for (i in names(markers_t)){
    print(i)
    temp <- markers_t[[i]][,4:ncol(markers_t[[i]])] %>% as.matrix()
    count <-  apply(temp, 1, function(x) sum(x > 1, na.rm = TRUE))
	pval = markers_t[[i]][,1]
	FDR = markers_t[[i]][,2]
    med_logFC <- apply(temp, 1, function(x) median(x, na.rm = TRUE))
    mean_logFC <- apply(temp, 1, function(x) mean(x, na.rm = TRUE))
    gene <- row.names(temp)
    cluster <- i
    markers_summary_t[[i]] <- cbind(gene, pval, FDR, count, med_logFC, mean_logFC, cluster) %>% as_tibble()
}
markers_summary_w <- markers_summary %>% bind_rows() %>% mutate(pval = as.numeric(pval), FDR = as.numeric(FDR), count = as.numeric(count), med_auc = as.numeric(med_auc), mean_auc = as.numeric(mean_auc))
markers_summary_t <- markers_summary_t %>% bind_rows() %>% mutate(pval = as.numeric(pval), FDR = as.numeric(FDR), count = as.numeric(count), med_logFC  = as.numeric(med_logFC), mean_logFC = as.numeric(mean_logFC))

save(markers_wilcox, markers_t, file = args[6])
save(markers_summary_w, markers_summary_t, file = args[7])
