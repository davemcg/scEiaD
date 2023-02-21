args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(SingleCellExperiment)
library(DESeq2)

load(args[1])
against <- str_extract(args[1], 'cluster|Cell\\w+')
meta <- colData(org_mat)
counts <- assay(org_mat)
colnames(counts) <- paste(meta$study_accession,
							meta[,against],
								sep = '__')

write_tsv(as_tibble(meta), file = gsub('diff_testing', 'site', args[1]) %>% gsub('deseq2obj.Rdata', 'meta.tsv.gz', .))
write_tsv(as_tibble(counts, rownames = 'Gene') , file = gsub('diff_testing', 'site', args[1]) %>% gsub('deseq2obj.Rdata', 'pseudoCounts.tsv.gz', .))
							

