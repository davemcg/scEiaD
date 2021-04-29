library(tidyverse)
human_t2g <- read_tsv('v35.tr2gX.tsv', col_names = F) 
mouse_t2g <- read_tsv('vM25.tr2gX.tsv',col_names =c('extra', 'transcript_id', 'gene_id' )) 
# this file also in the scEiaD git "data" dir
# so you can get it via
# wget https://github.com/davemcg/scEiaD/raw/master/data/ensembl_biomart_human2mouse_macaque.tsv
converter <- read_tsv('/data/mcgaugheyd/projects/nei/mcgaughey/scEiaD_me/references/ensembl_biomart_human2mouse_macaque.tsv') %>% 
  select(human_gid = `Gene stable ID`, gene_id = `Mouse gene stable ID`) %>% distinct %>% 
  filter(!duplicated(human_gid), !duplicated(gene_id)) %>% distinct
mouse_t2g_w_human <- left_join(mouse_t2g, converter) %>% 
  mutate(human_gid_w_mouse = replace(human_gid,is.na(human_gid), gene_id[is.na(human_gid)] ))
mouse_t2g_w_human %>% 
	select(extra, transcript_id, human_gid_w_mouse) %>% 
	write_tsv('vM25.tr2gX.humanized.tsv', col_names = F)
