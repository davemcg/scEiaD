library(tidyverse)
# quantify_performance.R

# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
# k-nearest neighbor batch-effect test (kBET) [22], local inverse Simpson’s index 
# (LISI) [13, 23], average silhouette width (ASW) [24], and adjusted rand index 
# (ARI) benchmarking metrics [25].

args <- commandArgs(trailingOnly = TRUE)
# umap
load(args[1])
# full obj
load(args[2])

# find proper embedding to use
reductions <- names(integrated_obj@reductions)
reductions <- reductions[!grepl('UMA', reductions, ignore.case = TRUE)]
if (length(reductions) == 1 && reductions == 'PCA'){
	reduction <- 'PCA'
} else { reduction = reductions[reductions != 'PCA'][1]
}
scores <- list()
# kBET
# devtools::install_github('theislab/kBET')
# https://github.com/theislab/kBET
set.seed(12534)
# take up to 3000 from each organism/celltype
cutdown <- umap %>% 
  rowid_to_column('ID') %>% 
  filter(!CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !is.na(CellType)) %>% 
  group_by(organism, CellType) %>% 
  sample_n(3000, replace = TRUE) %>% 
  unique()
silhouette <- function(obj){
  indices <- obj$ID
  faux_pca <- list()
  faux_pca$x <- integrated_obj@reductions[[reduction]]@cell.embeddings[indices,]
  batch <- umap[indices, 'batch'] %>% pull(1) %>% as.factor() %>% as.numeric()
  dims = ncol(integrated_obj@reductions[[reduction]]@cell.embeddings[indices,])
  kBET::batch_sil(faux_pca, batch, nPCs = dims)
}
# have to subsample as this requires making a dist obj! 
# -1 overlapping (good if you are evaluating batch mixing)
# 1 distinct clusters (bad if you are evaluating batch)
# 0 is random
scores$silhouette <- silhouette(cutdown)

# run silhouette by cell type
scores$silhouette_CellType <- cutdown %>% 
  group_by(CellType) %>% group_map(~ silhouette(.x))

# batch <- umap[cutdown,'batch'] %>% pull(1) %>% as.factor() %>% as.numeric()
# not much point running kBET for me as it returns fail for everything
# batch.estimate <- kBET::kBET(faux_pca$x, batch)


# LISI
# https://github.com/immunogenomics/LISI
# devtools::install_github("immunogenomics/lisi")
# generates a score PER CELL which is a metric of number of types of neighbors
# lower is better (pure cell population within region)
scores$LISI <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdown$ID,], 
                                  cutdown,
                                  'CellType')

# ARI
# https://davetang.org/muse/2017/09/21/adjusted-rand-index/
# very slow
# install.packages('clues')
ari <- function(obj){
  clues::adjustedRand(obj$cluster, obj$CellType %>% as.factor() %>% as.numeric())
}
scores$RI <- ari(cutdown)
# by species
scores$RI_species <- cutdown %>% 
  group_by(organism) %>% 
  group_map(~ ari(.x))

scores$Barcode <- cutdown %>% ungroup() %>% select(Barcode)
scores$umap_cutdown <- cutdown

save(scores, file = args[3])
