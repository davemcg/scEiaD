# extract UMAP coords labelled with meta data
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(tidyverse)
# load seurat obj
load(args[1])

# give seurat obj name as input
# assign to `obj`
obj <- get(args[2])

# load meta data (including inferred cell types)
load(args[3])
orig_meta <- obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(obj[['umap']]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., meta %>% select(Barcode, Age:new_CellType_transfer))

saveRDS(umap, file = args[4])


