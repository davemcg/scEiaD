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
# load cluster
load(args[3])
umap$cluster <- meta %>% pull(2)

print(args[2])

if (grepl('onlyWELL', args[2])){
	cutdown <- umap %>% rowid_to_column('ID')
	cutdownSUB <- umap %>% rowid_to_column('ID')
	scores <- list()
	reductions <- names(integrated_obj@reductions)
    reductions <- reductions[!grepl('UMA', reductions, ignore.case = TRUE)]
    if (length(reductions) == 1 && reductions == 'PCA'){
        reduction <- 'PCA'
    } else { reduction = reductions[reductions != 'PCA'][1]
    }
} else {
	umap <- umap %>%  mutate(SubCellType = gsub('p_','', SubCellType)) %>% mutate(SubCellType = gsub('f_','', SubCellType))
	umap$SubCellType[grepl('^RB|Rods|Peri|^MG$|^Mic$', umap$SubCellType)] <- NA

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
	# take up to 2000 from each organism/celltype
	print('cutdown begin')


	cutdown <- umap %>% 
	  rowid_to_column('ID') %>% 
	  filter(!CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
	         !is.na(CellType)) %>% 
	  group_by(organism, CellType) %>% 
	  sample_n(2000, replace = TRUE) %>% 
	  unique()
	# remove celltypes which have fewer than 100 cells
	keep <- umap %>% group_by(CellType) %>% summarise(count = n()) %>% filter(count > 99) %>% pull(CellType)
	cutdown <- cutdown %>% filter(CellType %in% keep)
	# cutdown against subcelltype
	cutdownSUB <- umap %>% 
	  rowid_to_column('ID') %>% 
	  filter(!is.na(SubCellType)) %>% 
	  group_by(organism, SubCellType) %>% 
	  sample_n(2000, replace = TRUE) %>% 
	  unique()
	# remove celltypes which have fewer than 100 cells
	keep <- umap %>% group_by(SubCellType) %>% summarise(count = n()) %>% filter(count > 99) %>% pull(SubCellType)
	cutdownSUB <- cutdownSUB %>% filter(SubCellType %in% keep)
}
# kBET silhouette function
# user can select what silhoette is calcualted AGAINST 
# have to subsample input as this requires making a dist obj! 
# -1 overlapping (good if you are evaluating batch mixing)
# 1 distinct clusters (bad if you are evaluating batch)
# 0 is random
silhouette <- function(obj, against = 'batch'){
  out <- 0
  try({
    indices <- obj$ID
    faux_pca <- list()
    faux_pca$x <- integrated_obj@reductions[[reduction]]@cell.embeddings[indices,]
    
	metric <- umap[indices, against] %>% pull(1) %>% as.factor() %>% as.numeric()
    dims = ncol(integrated_obj@reductions[[reduction]]@cell.embeddings[indices,])
    out <- kBET::batch_sil(faux_pca, metric, nPCs = dims)
  })
  out
}
print('scoring starts')
print('silhouette')
scores$silhouette_batch <- silhouette(cutdown)
if (grepl('onlyDROPLET', args[2])){
	scores$silhouette_celltype <- silhouette(cutdown, against = 'CellType')
	scores$silhouette_subcelltype <- silhouette(cutdownSUB, against = 'SubCellType')
}

scores$silhouette_cluster <- silhouette(cutdown, against = 'cluster')
# run silhouette by cell type
#print('silhouette by celltype')
#scores$silhouette_CellType <- cutdown %>% 
#  group_by(CellType) %>% group_map(~ silhouette(.x))

# batch <- umap[cutdown,'batch'] %>% pull(1) %>% as.factor() %>% as.numeric()
# not much point running kBET for me as it returns fail for everything
# batch.estimate <- kBET::kBET(faux_pca$x, batch)


# LISI
# https://github.com/immunogenomics/LISI
# devtools::install_github("immunogenomics/lisi")
# generates a score PER CELL which is a metric of number of types of neighbors
# lower is better (pure cell population within region)
print('lisi')
if (grepl('onlyDROPLET', args[2])){	
	try({scores$LISI_subcelltype <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdownSUB$ID,], 
                                  cutdownSUB,
                                  'SubCellType') })
	scores$LISI_celltype <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdown$ID,], 
                                  cutdown,
                                  'CellType')
}
scores$LISI_batch <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdown$ID,],
                                  cutdown,
                                  'batch')
scores$LISI_cluster <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdown$ID,],
                                  cutdown,
                                  'cluster')
# ARI
# https://davetang.org/muse/2017/09/21/adjusted-rand-index/
# very slow
# install.packages('clues')
#print('ari')
#ari <- function(obj){
#  out <- 'fail'
#  try({
#    clues::adjustedRand(obj$cluster %>% as.character() %>% as.numeric(), obj$CellType %>% as.factor() %>% as.numeric())
#  })
#}
#scores$RI <- ari(cutdown)
#print('scoring over')
# by species
#scores$RI_species <- cutdown %>% 
#  group_by(organism) %>% 
#  group_map(~ ari(.x))

scores$Barcode <- cutdown %>% ungroup() %>% select(Barcode)
scores$umap_cutdown <- cutdown

save(scores, file = args[4])
