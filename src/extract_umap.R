# extract UMAP coords labelled with meta data
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(tidyverse)
# load integrated seurat obj
load(args[1])
# load labelled cell data
load(args[2])
# load predited cell data (+ labelled)
load(args[3])
# method
method <- args[5]

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
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'combatUMAP_'
} else {
  print(paste0("Why did you pick ", method, "?"))
}
reduction.name <- gsub('_','', reduction.key)


# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(integrated_obj[[reduction.name]]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., cell_info_labels %>% dplyr::rename(Barcode = value)) %>% 
  left_join(., predictions %>% as_tibble(rownames = 'Barcode') %>% select(Barcode, CellType_predict = `predicted.id`))

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
# Vsx1, Grm6, Cabp5 for bipolar cells
# Elavl4 for RGC
# CX3CR1 for microglia
core_markers <- c('Rho','Opn1sw', 'Sfrp2', 'Hes6', 'Tfap2a', 'Isl1', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4', 'Best1', 'Grm6', 'Cabp5', "Cx3cr1")
core_expression <- FetchData(integrated_obj, toupper(core_markers)) %>% as_tibble(rownames = 'Barcode')

umap <- left_join(umap, avg_marker, by = 'Barcode') %>% 
  left_join(., core_expression, by = 'Barcode')

save(umap, file = args[4])


