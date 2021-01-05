library(tidyverse)
library(pool)
library(RSQLite)

args = commandArgs(trailingOnly=TRUE)


scEiaD <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)

load(args[2])
load(args[3])

PB_results <- bind_rows(PB_resultsC2, PB_resultsABC) %>% arrange(FDR, -abs(logFC))
# update gene name
gene_id_converter <- read_tsv('references/ensembl_biomart_human2mouse_macaque.tsv', skip = 1,
                              col_names= c('hs_gene_id','hs_gene_id_v', 'mm_gene_id', 'mf_gene_id',
                                           'hs_gene_name', 'mf_gene_name', 'mm_gene_name')) %>%
  select(-hs_gene_id_v)
PB_results <- PB_results %>% left_join(gene_id_converter %>% select(hs_gene_id, hs_gene_name) %>% unique(), by = c('Gene' = 'hs_gene_id')) %>% mutate(Gene = paste0(hs_gene_name, ' (', Gene, ')'))  %>% select(-hs_gene_name)

dbWriteTable(scEiaD, 'PB_results', PB_results, overwrite = TRUE)
db_create_index(scEiaD, table = 'PB_results', columns = c('Gene'))
db_create_index(scEiaD, table = 'PB_results', columns = c('test'))
db_create_index(scEiaD, table = 'PB_results', columns = c('PB_Test'))

# doublets
load(args[4])
dbWriteTable(scEiaD, 'doublets', doublet_call_table, overwrite = TRUE)
db_create_index(scEiaD, table = 'doublets', columns = c('Barcode'))

# extract pseudobulk tests and terms 
PB_Test_terms <- scEiaD %>% 
  tbl('PB_results') %>% 
  select(PB_Test) %>% 
  distinct() %>% 
  as_tibble()

terms <- list()
for(i in PB_Test_terms %>% pull(1)){
  print(i)
  terms[[i]] <- scEiaD %>% 
    tbl('PB_results') %>% 
    filter(PB_Test == i) %>% 
    select(test) %>% 
    distinct() %>% 
    pull(1) %>% sort() %>% 
    paste(collapse = '___')
}

PB_Test_terms$terms <- terms %>% unlist()
dbWriteTable(scEiaD, name = 'PB_Test_terms', PB_Test_terms, overwrite = TRUE)


# update meta table
meta_filter <- left_join(scEiaD %>% tbl('metadata_filter'),
                         scEiaD %>% tbl('doublets'), by ='Barcode') %>%
  collect() %>%
  mutate(`Doublet Probability` = as.numeric(`Doublet Probability`),
         doublet_score_scran = as.numeric(doublet_score_scran)) %>%
  mutate(PMID = as.character(PMID)) %>%
  mutate(TechType = case_when(Platform %in% c('10xv2','10xv3','DropSeq') ~ 'Droplet',
                          TRUE ~ 'Well')) 
# output for user download
write_tsv(meta_filter, path = 'site/metadata_filter.tsv.gz')

# coordinates
# get coords for cell labels
celltype_predict_labels <-  meta_filter %>%
              group_by(CellType_predict) %>%
              summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) 

celltype_labels <- meta_filter %>%
  group_by(CellType) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
# tabulamuris_labels <- meta_filter %>%
#   group_by(TabulaMurisCellType) %>%
#   summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
tabulamuris_predict_labels <- meta_filter %>%
  group_by(TabulaMurisCellType_predict) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
# get coords for cell labels
cluster_labels <- meta_filter %>%
              group_by(cluster) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# haystack
load(args[5])
haystack <- scH$results %>% as_tibble(rownames = 'hs_gene_id') %>% left_join(gene_id_converter %>% select(hs_gene_id, hs_gene_name) %>% unique()) %>% mutate(Gene = paste0(hs_gene_name, ' (', hs_gene_id, ')')) %>% select(Gene, D_KL, `log.p.vals`, `log.p.adj`, `T.counts`)


dbWriteTable(scEiaD, 'metadata_filter', meta_filter, overwrite = TRUE)
dbWriteTable(scEiaD, 'celltype_predict_labels', celltype_predict_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'celltype_labels', celltype_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'tabulamuris_predict_labels', tabulamuris_predict_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'cluster_labels', cluster_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'haystack', haystack, overwrite = TRUE)
db_create_index(scEiaD, table = 'haystack', columns = c('Gene'))
