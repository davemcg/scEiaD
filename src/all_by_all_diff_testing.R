library(Seurat)
library(tidyverse)
library(future)
library(furrr)

args <- commandArgs(trailingOnly = TRUE)

plan("multiprocess", workers = as.numeric(args[1]))
options(future.globals.maxSize = 500000 * 1024^2)




# cluster
load(args[2])
# integrated_obj
load(args[3])
# output obj
out <- args[4]

Idents(integrated_obj) <- meta[,2] %>% pull(1)

cluster_names <- meta[,2] %>%
  pull(1) %>%
  unique()

calc_diff_combn <- function(i, obj, cluster_IDs){
  first <- combn(cluster_IDs,2)[1,i]
  second <- combn(cluster_IDs,2)[2,i]
  comparison = paste0(first,'__',second)
  print(comparison)
  out <- FindMarkers(obj, ident.1 = first, ident.2 = second, logfc.threshold = 0.5, verbose = FALSE)
  out$comparison <- comparison
  out
}

calc_diff <- function(i, obj, cluster_IDs){
  out <- FindMarkers(obj, ident.1 = cluster_IDs[i], logfc.threshold = 0.5, verbose = FALSE)
  out$comparison <- cluster_IDs[i] 
  if (nrow(out) > 0) {
  	out
  } else {'no significant results'}
}

de_rna_combn <- future_map(1:ncol(combn(cluster_names,2)), 
					calc_diff_combn, 
					integrated_obj, 
					cluster_names, 
					.progress = TRUE)
de_rna_all <- future_map(1:length(cluster_names), 
				calc_diff, 
				integrated_obj, 
				cluster_names, 
				.progress = TRUE)

save(de_rna_combn, de_rna_all, file = out_file)

