Sys.setenv(RETICULATE_PYTHON = "/data/mcgaugheyd/conda/bin/python3.7")
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

save(doublet_call_table, file = args[2])
