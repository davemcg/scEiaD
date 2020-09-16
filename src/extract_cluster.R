library(Seurat)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])

meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% select(Barcode, contains('cluster'))
graph <- integrated_obj@misc[[names(integrated_obj@misc)]]
# swap scArches subcluster for cluster as default params get SO FEW DAMN CLUSTERS 
if (grepl('scArches', args[1])){
	meta <- meta[,c(1,3,2)]
}
save(meta, file = args[2])
save(graph, file = args[3])
