library(tidyverse)
library(dtplyr)
human_mouse <- read_tsv('~/git/scEiaD/data/human_mouse_hcop_fifteen_column.txt.gz')
hmm <- human_mouse %>% select(human_name, human_ensembl_gene, human_symbol, mouse_ensembl_gene, mouse_symbol) %>% distinct()

human_chicken <- read_tsv('~/git/scEiaD/data/human_chicken_hcop_fifteen_column.txt.gz')
hch <- human_chicken %>% select(human_name, human_ensembl_gene, human_symbol, chicken_ensembl_gene, chicken_symbol)

human_zebrafish <- read_tsv('~/git/scEiaD/data/human_zebrafish_hcop_fifteen_column.txt.gz')
hzf <- human_zebrafish %>% select(human_name, human_ensembl_gene, human_symbol, zebrafish_ensembl_gene, zebrafish_symbol)

human_macaque <- read_tsv('~/git/scEiaD/data/human_macaque_hcop_fifteen_column.txt.gz')
hmf <- human_macaque %>% select(human_name, human_ensembl_gene, human_symbol, macaque_ensembl_gene, macaque_symbol)


#human_fly <- read_tsv('~/Downloads/human_fruitfly_hcop_fifteen_column.txt.gz')
#hfly <- human_fly %>% select(human_name, human_ensembl_gene, human_symbol, `fruit fly_ensembl_gene`, `fruit fly_symbol`)

gene_id_converter <- lazy_dt(hmm) %>% 
  left_join(lazy_dt(hch), by = c('human_name','human_ensembl_gene','human_symbol')) %>% 
  left_join(lazy_dt(hzf), by = c('human_name','human_ensembl_gene','human_symbol')) %>% 
  left_join(lazy_dt(hmf), by = c('human_name','human_ensembl_gene','human_symbol'))

gene_id_converter <- as_tibble(gene_id_converter)


colnames(gene_id_converter) <- c('human_name','hs_gene_id',
                                 'hs_gene_name',
                                 'mf_gene_id', 'mf_gene_name',
                                 'zf_gene_id', 'zf_gene_name',
                                 'gg_gene_id', 'gg_gene_name',
                                 'mm_gene_id', 'mm_gene_name')
