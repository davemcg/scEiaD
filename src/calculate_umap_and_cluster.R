library(tidyverse)
library(Seurat)
library(future)
plan(strategy = "multicore", workers = 4)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
max_dims = args[2] %>% as.numeric()
# cell labels
load(args[3])
# integrated_obj
load(args[4])


create_umap_and_cluster <- function(integrated_obj, 
                                    max_dims = 20, 
                                    reduction = 'pca',
                                    reduction.name = 'ccaUMAP',
                                    reduction.key = 'ccaUMAP_'){
  # UMAP
  integrated_obj <- RunUMAP(integrated_obj, 
                            dims = 1:max_dims, 
                            reduction = reduction, 
                            reduction.name = reduction.name,
                            reduction.key = reduction.key)
  # clustering 
  integrated_obj <- FindNeighbors(integrated_obj, 
                                  reduction = reduction,
                                  dims = 1:max_dims, 
                                  nn.eps = 0.5)
  integrated_obj <- FindClusters(integrated_obj, 
                                 resolution = c(0.1,0.3,0.6,0.8,1,2,3,4,5),
                                 save.SNN = TRUE,
                                 do.sparse = TRUE,
                                 algorithm = 2,
                                 random.seed = 23)
  integrated_obj
}

if (method == 'CCA'){
  reduction <- 'pca'
  reduction.key <- 'ccaUMAP_'
} else if (method == 'scanorama'){
  reduction <- 'scanorama'
  reduction.key <- 'scanoramaUMAP_'
} else if (method == 'harmony'){
  reduction <- 'harmony'
  reduction.key <- 'harmonyUMAP_'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
  reduction.key <- 'mnnUMAP_'
} else if (method == 'none'){
  reduction <- 'pca'
  reduction.key <- 'noneUMAP_'
} else if (method == 'combat'){
  reduction <- 'pca'
  reduction.key <- 'combatUMAP_'
} else {
  print(paste0("Why did you pick ", method, "?"))
}
reduction.name <- gsub('_','', reduction.key)
integrated_obj <- create_umap_and_cluster(integrated_obj, 
                                          max_dims,
                                          reduction,
                                          reduction.name = reduction.name,
                                          reduction.key = reduction.key)

save(integrated_obj, file = args[5], compress = FALSE )

# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(integrated_obj[[reduction.name]]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., cell_info_labels %>% dplyr::rename(Barcode = value))
umap$Method <- method
colnames(umap)[2:3] <- c('UMAP_1', 'UMAP_2')

# pull eye markers
msig <- readxl::read_xls('~/git/massive_integrated_eye_scRNA/data/scsig.v1.0.metadata.xls')
msig <- msig %>% filter(grepl('Visual', `Source Organ System`)) %>% select(`Gene Set Standard Name`, `Source Organism`, `Exact Source In Publication`, `PubMed ID`, `Additional Details`, `Raw Member IDs from Source Publication`) %>% mutate(Gene = str_split(`Raw Member IDs from Source Publication`, ',')) %>% unnest() 
msig <- msig %>% mutate(Num = row.names(msig))
# filter to genes that are in the scale.data slot in the DefaultAssay 
assay <- DefaultAssay(integrated_obj)
msig_present <- msig %>% filter(Gene %in% row.names(integrated_obj@assays[[assay]]@scale.data)) %>% 
  filter(!grepl('MT', Gene)) # no mito genes
top <- msig_present %>% group_by(`Additional Details`) %>% top_n(., 10, wt = Num)
# hybrid expresions system where I average the top 10 genes in each marker group
expression <- list()
for (i in top$`Additional Details` %>% unique()){
  expression[[i]] <- FetchData(integrated_obj, top %>% filter(`Additional Details` == i) %>% pull(Gene))
}
avg_marker <- expression %>% map(rowMeans) %>% map(enframe) %>% bind_rows(.id = 'id')
colnames(avg_marker) <- c('CellType', 'Barcode', 'ScaleData')
avg_marker <- avg_marker %>% spread(CellType, ScaleData)

# add core markers
# Rho for rods
# Opn1sw for cones
# Sfrp2 for early progenitors
# Hes6 for neurogenic
# Tfap2a for amacrine
# Ccnd1 for late progenitors
# Aqp4 for muller glia
# Vsx1 to bipolar
# Elavl4 for RGC
core_markers <- c('Rho','Opn1sw', 'Sfrp2', 'Hes6', 'Tfap2a', 'Isl1', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4', 'Best1')
core_expression <- FetchData(integrated_obj, toupper(core_markers)) %>% as_tibble(rownames = 'Barcode')

umap <- left_join(umap, avg_marker, by = 'Barcode') %>% 
  left_join(., core_expression, by = 'Barcode')

save(umap, file = args[6])
