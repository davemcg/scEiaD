# quantify_performance.R

# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
# k-nearest neighbor batch-effect test (kBET) [22], local inverse Simpson’s index 
# (LISI) [13, 23], average silhouette width (ASW) [24], and adjusted rand index 
# (ARI) benchmarking metrics [25].

# kBET
# https://github.com/theislab/kBET

# ASW
# https://github.com/theislab/kBET#compute-a-silhouette-width-and-pca-based-measure

# LISI
# https://github.com/immunogenomics/LISI
# devtools::install_github("immunogenomics/lisi")
lisi_res <- lisi::compute_lisi(integrated_obj@reductions$scVI@cell.embeddings, 
                               umap$CellType_predict)

# ARI
# https://davetang.org/muse/2017/09/21/adjusted-rand-index/
# very slow
# clues::adjustedRand(umap$cluster, umap$CellType_predict %>% as.factor() %>% as.numeric())