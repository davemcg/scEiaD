library(tidyverse)
library(glue)

source(glue('{git_dir}src/make_gene_id_converter_table.R'))

maker <- function(cols){
	out <- gene_id_converter %>%  select(contains(cols)) %>% unique()
	out
}

hs_mm <- maker(c("hs_gene_id","hs_gene_name", "mm_gene_id", "mm_gene_name"))
hs_mf <- maker(c("hs_gene_id","hs_gene_name", "mf_gene_id", "mf_gene_name"))
hs_gg <- maker(c("hs_gene_id","hs_gene_name", "gg_gene_id", "gg_gene_name"))

write_tsv(hs_mm, 'site/hs_mm_converter.tsv.gz')
write_tsv(hs_mf, 'site/hs_mf_converter.tsv.gz')
write_tsv(hs_gg, 'site/hs_gg_converter.tsv.gz')
