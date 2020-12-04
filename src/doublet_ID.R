conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
library(glue)
Sys.setenv(RETICULATE_PYTHON = glue("{conda_dir}/bin/python3.7"))

scr <- reticulate::import('scrublet')
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
load(args[1])

doublet_calls <- list()
for (each_batch in unique(integrated_obj@meta.data$batch)){
	print(each_batch)	
	bc <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% filter(batch == each_batch) %>% pull(Barcode)
	matrix <- integrated_obj@assays$RNA@counts[,bc] %>% as.matrix()
	prep <- scr$Scrublet(matrix %>% t())
	out <- list()
	out <- tryCatch(prep$scrub_doublets(), error = function(e) {list(0,'FAIL')})
	out[[3]] <- bc
	df <- cbind(out[[3]], out[[1]], out[[2]])
	colnames(df) <- c('Barcode', 'Doublet Probability', 'Doublet')
	doublet_calls[[each_batch]] <- df
	
}
doublet_call_table <- doublet_calls %>% purrr::map(as_tibble) %>% bind_rows()

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(tidyverse)

int_sce <-  as.SingleCellExperiment(integrated_obj)
batches <- unique(int_sce$batch)
doublets <- list()
for (i in batches){
  print(i)
  doublet_score_scran <- scran::doubletCells(int_sce[,int_sce$batch == i])
  Barcode <- int_sce[,int_sce$batch == i] %>% colnames()
  out <- cbind(Barcode, doublet_score_scran) %>% as_tibble()
  doublets[[i]] <- out
}
doubl <- bind_rows(doublets)

doublet_call_table <- left_join(doublet_call_table, doubl)


save(doublet_call_table, file = args[2])
