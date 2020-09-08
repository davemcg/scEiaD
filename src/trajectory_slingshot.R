library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
load(args[1]) #load('/Volumes/OGVFB_BG/scEiaD/seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umap.Rdata')
load(args[2]) # load('/Volumes/OGVFB_BG/scEiaD/cluster/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__knn7.cluster.Rdata') 
load(args[3]) # load('/Volumes/OGVFB_BG/scEiaD/cell_info_labelled.Rdata')
load(args[4]) # load('/Volumes/OGVFB_BG/scEiaD/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umapFilter.predictions.Rdata')
start_clus <- as.integer(args[5])
num_cells <- as.integer(args[6])
reduction <- args[7]
out <- args[8]

integrated_obj@meta.data$cluster <- meta[,2] %>% pull(1)

orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')

umap <- Embeddings(integrated_obj[[reduction]]) %>% as_tibble(rownames = 'Barcode')  %>% 
  left_join(., integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% select(Barcode, cluster))


# Got: 
# Error in solve.default(s1 + s2) : system is computationally singular: reciprocal condition number = 1.16521e-18
# trying to add some noise as per: https://github.com/kstreet13/slingshot/issues/35
#sce <- Seurat::as.SingleCellExperiment(integrated_obj)
sce <- SingleCellExperiment(list(counts = integrated_obj@assays$RNA@counts[,umap$Barcode]),
                                 colData= umap %>% select(cluster))


reducedDims(sce) <- list(scVI= umap %>% select(contains(reduction)) )


set.seed(2534)
sds_in <- apply(reducedDim(sce, "scVI"), 2, function(x) jitter(x))

sds <- slingshot(sds_in, 
                 clusterLabels = sce$cluster, 
                 start.clus = start_clus, 
                 stretch = 0, approx_points = 10000)
save(sds, sds_in, file = args[6], compress = FALSE)
