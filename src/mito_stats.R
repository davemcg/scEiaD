args= commandArgs(trailingOnly = T)
library(tidyverse)
wd = args[1]
drop_mt_pfx = args[2]
quant_path =args[3]
mito_gene_stem = args[4]
outfile = args[5]


files_drop <- list.files('pipeline_data/', pattern = drop_mt_pfx, recursive=TRUE, full.names=TRUE)
files_well <- list.files(quant_path, pattern = 'abundance_spliced.tsv.gz', recursive=TRUE, full.names=TRUE)
mito_gene_files <- list.files('references/', pattern = mito_gene_stem, full.names = T, recursive = T)
drop_mito <- lapply(files_drop, read_tsv) %>% bind_rows
mito_genes <- lapply(mito_gene_files, function(x) scan(x, character(), sep= '\n')) %>%  reduce( c)



well_meta <- list()
inc <- 1
for (i in files_well){
	print(inc); inc = inc + 1
	sample <- str_extract(i, 'SRS\\d+')
    mat <- read_tsv(i)
	m <- mat[,3, drop = FALSE] %>% as.matrix()
	row.names(m) <- mat$target_id
	colnames(m) <- sample
    mito_sum <- m[mito_genes[mito_genes %in% row.names(m)],1] %>% sum()
	all_sum <- m[,1] %>% sum()
	perc <- (mito_sum/all_sum) * 100
	well_meta[[sample]] <- perc

}



well_mito <- well_meta %>% bind_rows() %>% pivot_longer(everything()) %>% arrange(-value)
colnames(well_mito) <- c('Barcode', 'percent.mt')

mito <- bind_rows(well_mito, drop_mito %>% select(Barcode = barcode,`percent.mt` ))

write_tsv(mito, file = 'mito_counts.tsv')
