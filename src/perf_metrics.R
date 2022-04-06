library(tidyverse)
# quantify_performance.R

# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9
# k-nearest neighbor batch-effect test (kBET) [22], local inverse Simpson’s index 
# (LISI) [13, 23], average silhouette width (ASW) [24], and adjusted rand index 
# (ARI) benchmarking metrics [25].

args <- commandArgs(trailingOnly = TRUE)
# save(args, file = 'testing/perf.Rdata')
# print('asd')
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
	umap$CellType <- gsub('Early ','',umap$CellType)
	umap$CellType <- gsub('Late  ','',umap$CellType)
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
	  filter(CellType != 'RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cell',
			 !grepl('Corneal', CellType),
			 !CellType %in% c('Doublet', 'Doublets'),
	         !is.na(CellType)) %>% 
	  group_by(study_accession, CellType) %>% 
	  sample_n(300, replace = TRUE) %>% 
	  unique()
	# remove celltypes which have fewer than 100 cells
	keep <- umap %>% group_by(CellType) %>% summarise(count = n()) %>% filter(count > 99) %>% pull(CellType)
	cutdown <- cutdown %>% filter(CellType %in% keep)
	# cutdown against subcelltype
	cutdownSUB <- umap %>% 
	  rowid_to_column('ID') %>% 
	  filter(!is.na(SubCellType)) %>% 
	  filter(!grepl('Cone|Rod|MG|Muller|RBC|Endo|Peri', SubCellType)) %>% 
      mutate(SubCellType = gsub('f_', '', SubCellType)) %>% mutate(SubCellType = gsub('p_', '', SubCellType)) %>%
	  group_by(organism, SubCellType) %>% 
	  sample_n(1000, replace = TRUE) %>% 
	  unique()
	# remove celltypes which have fewer than 100 cells
	keep <- umap %>% group_by(SubCellType) %>% summarise(count = n()) %>% filter(count > 99) %>% pull(SubCellType)
	cutdownSUB <- cutdownSUB %>% filter(SubCellType %in% keep)
	
	#cutdown BIG
	cutdownBIG <- umap %>% 
	  rowid_to_column('ID') %>% 
	  filter(CellType != 'RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cell',
			 !grepl('Corneal', CellType),
			 !CellType %in% c('Doublet', 'Doublets'),
	         !is.na(CellType)) %>% 
	  group_by(study_accession, CellType) %>% 
	  sample_n(10000, replace = TRUE) %>% 
	  unique()
	# remove celltypes which have fewer than 100 cells
	keep <- umap %>% group_by(CellType) %>% summarise(count = n()) %>% filter(count > 99) %>% pull(CellType)
	cutdownBIG <- cutdownBIG %>% filter(CellType %in% keep)
	
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
cluster_count <- meta %>% pull(2) %>% table() %>% as_tibble()
colnames(cluster_count) <- c('Cluster','Count')
scores$cluster_count <- cluster_count 
print('scoring starts')
print('silhouette')
print(dim(cutdown))
scores$silhouette_batch <- silhouette(cutdown)
#if (!grepl('onlyWELL', args[2])){
#	scores$silhouette_celltype <- silhouette(cutdown, against = 'CellType')
#	scores$silhouette_subcelltype <- silhouette(cutdownSUB, against = 'SubCellType')
#}

scores$silhouette_cluster <- silhouette(cutdown, against = 'cluster')
# run silhouette by cell type
print('silhouette by celltype by celltype')
scores$silhouette_CellType_groupBy <- cutdownBIG %>% 
  group_by(CellType) %>% group_map(~ silhouette(.x))
names(scores$silhouette_CellType_groupBy) <- cutdownBIG %>% group_by(CellType) %>% summarise(Count = n()) %>% pull(CellType)
# batch <- umap[cutdown,'batch'] %>% pull(1) %>% as.factor() %>% as.numeric()
# not much point running kBET for me as it returns fail for everything
# batch.estimate <- kBET::kBET(faux_pca$x, batch)


# LISI
# https://github.com/immunogenomics/LISI
# devtools::install_github("immunogenomics/lisi")
# generates a score PER CELL which is a metric of number of types of neighbors
# lower is better (pure cell population within region)
print('lisi')
LISI_info <- function(df, cutdown_obj){
	df$Barcode <- cutdown_obj$Barcode
	df
}
if (!grepl('onlyWELL', args[2])){	
	scores$LISI_subcelltype <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdownSUB$ID,], 
                                  cutdownSUB,
                                  'SubCellType')
	scores$LISI_celltype <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdownBIG$ID,], 
                                  cutdownBIG,
                                  'CellType')
	scores$LISI_subcelltype <- LISI_info(scores$LISI_subcelltype, cutdownSUB)
	scores$LISI_celltype <- LISI_info(scores$LISI_celltype, cutdownBIG)
}
scores$LISI_batch <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdownBIG$ID,],
                                  cutdownBIG,
                                  'batch')
scores$LISI_cluster <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdownBIG$ID,],
                                  cutdownBIG,
                                  'cluster')

scores$LISI_batch <- LISI_info(scores$LISI_batch, cutdownBIG)
scores$LISI_cluster <- LISI_info(scores$LISI_cluster, cutdownBIG)

#scores$LISI_celltype <- lisi::compute_lisi(integrated_obj@reductions[[reduction]]@cell.embeddings[cutdown$ID,],
#                                  cutdown,
#                                  'CellType')
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
scores$Barcode_big <- cutdownBIG %>% ungroup() %>% select(Barcode)
scores$umap_big <- cutdownBIG
save(scores, file = args[4])
