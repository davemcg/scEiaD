library(Seurat)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])

meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% select(Barcode, contains('cluster'))
graph <- integrated_obj@misc[[names(integrated_obj@misc)]]
save(meta, file = args[2])
save(graph, file = args[3])
