library(furrr)
library(tictoc)
library(Seurat)
plan(strategy = "multicore", workers = 24)
options(future.globals.maxSize = 500000 * 1024^2)
args <- commandArgs(trailingOnly = TRUE)

# cluster 
load(args[1])
# integrated_obj
load(args[2])
# output obj
out <- args[3]

calc_diff <- function(i, obj){
  first <- combn(cluster_nums,2)[1,i]
  second <- combn(cluster_nums,2)[2,i]
  comparison = paste0(first,'__',second)
  print(comparison)
  out <- FindMarkers(obj, ident.1 = first, ident.2 = second, logfc.threshold = 0.5)
  out$comparison <- comparison
  out
}

cluster_nums <- meta[,2] %>%
  pull(1) %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>% 
  as.integer()

tic()
de_rna <- map(1:ncol(combn(cluster_nums,2)), calc_diff, integrated_obj)
toc()

save(de_rna,  file = out)

