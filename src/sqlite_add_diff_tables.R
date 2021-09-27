library(tidyverse)
library(pool)
library(RSQLite)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)
print(args)

git_dir = Sys.getenv('SCIAD_GIT_DIR')
scEiaD <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)


load(args[2]) # CellType_predict
diff_testing_CTp <- diff_testing


load(args[3]) # CellType
diff_testing_CT <- diff_testing

load(args[4]) # cluster
diff_testing_cluster <- diff_testing

wilcox_diff_testing <- bind_rows(diff_testing_CTp, diff_testing_CT, diff_testing_cluster) %>% arrange(-AUC)



# update gene name
print('update gene name')
#gene_id_converter <- read_tsv('references/ensembl_biomart_human2mouse_macaque.tsv', skip = 1,
#                              col_names= c('hs_gene_id','hs_gene_id_v', 'mm_gene_id', 'mf_gene_id',
#                                           'hs_gene_name', 'mf_gene_name', 'mm_gene_name')) %>%
#  dplyr::select(-hs_gene_id_v)
library(glue)
source(glue('{git_dir}src/make_gene_id_converter_table.R'))
#gene_id_converter <- read_tsv(glue('{git_dir}/data/ensembl_biomart_human2mouse_macaque_chick_ZF.tsv.gz'), skip = 1,
#                              col_names= c('hs_gene_id','hs_gene_id_v',
#                                           'gg_gene_id', 'gg_gene_name',
#                                           'mf_gene_id', 'mf_gene_name',
#                                           'mm_gene_id', 'mm_gene_name',
#                                           'zf_gene_id', 'zf_gene_name',
#                                           'hs_gene_name')) %>%
#  select(-hs_gene_id_v)
print('diff testing tables')
wilcox_diff_testing <- wilcox_diff_testing %>% left_join(gene_id_converter %>% dplyr::select(hs_gene_id, hs_gene_name) %>% unique(), by = c('Gene' = 'hs_gene_id')) %>% mutate(Gene = paste0(hs_gene_name, ' (', Gene, ')'))  %>% dplyr::select(-hs_gene_name)

dbWriteTable(scEiaD, 'wilcox_diff_testing', wilcox_diff_testing, overwrite = TRUE)
db_create_index(scEiaD, table = 'wilcox_diff_testing', columns = c('Gene'))
db_create_index(scEiaD, table = 'wilcox_diff_testing', columns = c('Base'))
db_create_index(scEiaD, table = 'wilcox_diff_testing', columns = c('Group'))
wilcox_group_base_sets <- wilcox_diff_testing %>% dplyr::select(Group, Base)  %>% unique() %>% arrange(Group, Base)
dbWriteTable(scEiaD, 'wilcox_diff_testing_sets', wilcox_group_base_sets, overwrite = TRUE)
wilcox_diff_testing_genes <- scEiaD %>% tbl('wilcox_diff_testing') %>% pull(Gene) %>% unique()
wilcox_diff_testing_genes <- wilcox_diff_testing_genes %>% tibble::enframe(value = 'Gene') %>% dplyr::select(-name) %>% arrange()
dbWriteTable(scEiaD, 'wilcox_diff_testing_genes', wilcox_diff_testing_genes, overwrite = TRUE)

# doublets
print('write doublets')
load(args[5])
dbWriteTable(scEiaD, 'doublets', doublet_call_table, overwrite = TRUE)
db_create_index(scEiaD, table = 'doublets', columns = c('Barcode'))

# extract pseudobulk tests and terms 
# wilcox_diff_terms <- scEiaD %>% 
#  tbl('wilcox_diff_testing') %>% 
#  dplyr::select(Group) %>% 
#  distinct() %>% 
#  as_tibble()

# terms <- list()
# for(i in wilcox_diff_terms %>% pull(1)){
#  print(i)
#  terms[[i]] <- scEiaD %>% 
#    tbl('wilcox_diff_testing') %>% 
#    filter(Group == i) %>% 
#    dplyr::select(test) %>% 
#    distinct() %>% 
#    pull(1) %>% sort() %>% 
#    paste(collapse = '___')
#}

#wilcox_diff_terms$terms <- terms %>% unlist()
#dbWriteTable(scEiaD, name = 'wilcox_diff_terms', wilcox_diff_terms, overwrite = TRUE)


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
#load(args[6])
#haystack <- scH$results %>% as_tibble(rownames = 'hs_gene_id') %>% left_join(gene_id_converter %>% dplyr::select(hs_gene_id, hs_gene_name) %>% unique()) %>% mutate(Gene = paste0(hs_gene_name, ' (', hs_gene_id, ')')) %>% dplyr::select(Gene, D_KL, `log.p.vals`, `log.p.adj`, `T.counts`)


dbWriteTable(scEiaD, 'metadata_filter', meta_filter, overwrite = TRUE)
dbWriteTable(scEiaD, 'celltype_predict_labels', celltype_predict_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'celltype_labels', celltype_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'tabulamuris_predict_labels', tabulamuris_predict_labels, overwrite = TRUE)
dbWriteTable(scEiaD, 'cluster_labels', cluster_labels, overwrite = TRUE)
#dbWriteTable(scEiaD, 'haystack', haystack, overwrite = TRUE)
#db_create_index(scEiaD, table = 'haystack', columns = c('Gene'))
