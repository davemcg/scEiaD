library(tidyverse)

#args = commandArgs(trailingOnly=TRUE)

#gene_id_converter <- read_tsv('references/ensembl_biomart_human2mouse_macaque.tsv', skip = 1,
#                              col_names= c('hs_gene_id','hs_gene_id_v', 'mm_gene_id', 'mf_gene_id',
#                                           'hs_gene_name', 'mf_gene_name', 'mm_gene_name')) %>%
#  select(-hs_gene_id_v)
gene_id_converter <- read_tsv(glue('{git_dir}/data/ensembl_biomart_human2mouse_macaque_chick_ZF.tsv.gz'), skip = 1,
                              col_names= c('hs_gene_id','hs_gene_id_v',
                                           'gg_gene_id', 'gg_gene_name',
                                           'mf_gene_id', 'mf_gene_name',
                                           'mm_gene_id', 'mm_gene_name',
                                           'zf_gene_id', 'zf_gene_name',
                                           'hs_gene_name')) %>%
  select(-hs_gene_id_v)

gene_length <- function(ref){
    gtf <- rtracklayer::readGFF(ref)
	if (grepl('ENSMUSG', gtf$gene_id[1])){
		gene_size <- gtf %>% filter(type == 'transcript') %>% 
				select(gene_id, start, end) %>%
				mutate(l = abs(as.numeric(end)-as.numeric(start)), 
					   gene_id = gsub('\\.\\d+', '', gene_id)) %>%
				left_join(gene_id_converter, by = c('gene_id' = 'mm_gene_id')) %>% 
				group_by(hs_gene_id) %>%
				summarise(size = mean(l))
	} else {
		gene_size <- gtf %>% filter(type == 'transcript') %>% 
				select(gene_id, start, end) %>%
				mutate(l = abs(as.numeric(end)-as.numeric(start)), 
					   gene_id = gsub('\\.\\d+', '', gene_id)) %>%
				group_by(gene_id) %>%
				summarise(size = mean(l))
	}
	vec <- gene_size$size
	names(vec) <- gene_size[,1] %>% pull(1)
	vec
}


