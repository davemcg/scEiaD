args= commandArgs(trailingOnly = T)
#save(args, file = 'testing/mito_stargs.Rdata')
library(parallel)
library(tidyverse)
library(glue)
wd = args[1]
drop_mt_pfx = args[2]
quant_path =args[3]
mito_gene_stem = args[4]
patterns <- scan(args[5], what = character(), sep='\n') %>% paste0(collapse = '|')
outfile = args[6]


files_drop <- list.files('pipeline_data/', pattern = drop_mt_pfx, recursive=TRUE, full.names=TRUE)
files_well <- list.files(quant_path, pattern = 'abundance.tsv.gz', recursive=TRUE, full.names=TRUE)
mito_gene_files <- list.files('references/', pattern = mito_gene_stem, full.names = T, recursive = T)
drop_mito <- lapply(files_drop, read_tsv) %>% bind_rows
mito_genes <- lapply(mito_gene_files, function(x) scan(x, character(), sep= '\n')) %>%  reduce( c)

read_well_pt_mito <- function(file){
  sample <- str_extract(file, glue('({patterns})\\d+') )
  mat <- read_tsv(file, )
  m <- mat[,3, drop = FALSE] %>% as.matrix()
  row.names(m) <- mat$target_id
  colnames(m) <- sample
  mito_sum <- m[mito_genes[mito_genes %in% row.names(m)],1] %>% sum()
  all_sum <- m[,1] %>% sum()
  perc <- (mito_sum/all_sum) * 100
  return(tibble(srs=sample,barcode= sample, `percent.mt` = perc))
}

well_meta <- lapply(files_well, read_well_pt_mito)

well_mito <- well_meta %>% bind_rows()
mito <- bind_rows(well_mito, drop_mito)

write_tsv(mito, file = 'mito_counts.tsv')

