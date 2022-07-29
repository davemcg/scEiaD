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



diff <- data.table::fread(args[2])

join_by_diff_test <- function(diff_table, group, against = 'All'){
  diff_sub_hs<- diff_table %>% filter(Against %in% against, Group == group, grepl('ENSG', Gene))
  diff_sub_NOThs<- diff_table %>% filter(Against %in% against, Group == group, !grepl('ENSG', Gene))
  
  out_hs <- left_join(as_tibble(diff_sub_hs), 
                      gene_id_converter %>% select(hs_gene_id, hs_gene_name, human_name) %>% unique(), 
                      by = c('Gene' = 'hs_gene_id')) %>% 
    as_tibble()%>% 
    dplyr::rename(gene_name = hs_gene_name)
  
  out_NOThs <- left_join(as_tibble(diff_sub_NOThs), 
                         gene_id_converter %>% select(mf_gene_id, mf_gene_name, human_name) %>% unique(), 
                         by = c('Gene' = 'mf_gene_id')) %>% 
    as_tibble() %>% 
    dplyr::rename(gene_name = mf_gene_name)
  
  bind_rows(out_hs, out_NOThs)
}

types  <- diff$Against %>% unique()

ctp_all <- join_by_diff_test(diff, 'CellType_predict', 'All')
ctp_pairwise <- join_by_diff_test(diff, 'CellType_predict', types[types != 'All'])



ct_all <- join_by_diff_test(diff, 'CellType', 'All')
ct_pairwise <- join_by_diff_test(diff, 'CellType', types[types != 'All'])

cluster_all <- join_by_diff_test(diff, 'cluster', 'All')
cluster_pairwise <- join_by_diff_test(diff, 'cluster', types[types != 'All'])

diff_data <- bind_rows(ctp_all,
                       ctp_pairwise,
                       ct_all,
                       ct_pairwise,
                       cluster_all,
                       cluster_pairwise)


diff_data <- diff_data %>% 
  mutate(Gene = paste0(gene_name, ' (', Gene, ')')) %>% select(-gene_name) %>% dplyr::rename(`Gene Name` = human_name) %>% arrange(padj) %>% 
  mutate(Group = case_when(Group == 'cluster' ~ 'Cluster',
                           Group == 'CellType_predict' ~ 'CellType (Predict)',
                           TRUE ~ Group))

dbWriteTable(scEiaD, 'diff_testing', diff_data, overwrite = TRUE)

db_create_index(scEiaD, table = 'diff_testing', columns = c('Gene'))
db_create_index(scEiaD, table = 'diff_testing', columns = c('Base'))
db_create_index(scEiaD, table = 'diff_testing', columns = c('Group'))
db_create_index(scEiaD, table = 'diff_testing', columns = c('Against'))

group_base <- diff_data %>% select(Group, Base) %>% unique()
group_base_against <- diff_data %>% select(Group, Base, Against) %>% unique() %>% filter(Against != 'All')
diff_genes <- diff_data$Gene %>% unique() %>% enframe() %>% select(Gene = value)



dbWriteTable(scEiaD, 'diff_testing_sets', group_base_against, overwrite = TRUE)
dbWriteTable(scEiaD, 'diff_testing_genes', diff_genes, overwrite = TRUE)



# doublets
print('write doublets')
load(args[3])
dbWriteTable(scEiaD, 'doublets', doublet_call_table, overwrite = TRUE)
db_create_index(scEiaD, table = 'doublets', columns = c('Barcode'))
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

# gene name id table
gic <- gene_id_converter %>% filter(human_name != '-')
hs <- gic %>% select(human_name, hs_gene_id) %>% unique()
mm <- gic %>% select(human_name, mm_gene_id) %>% unique()
mf <- gic %>% select(human_name, mf_gene_id) %>% unique()
gg <- gic %>% select(human_name, gg_gene_id) %>% unique()

name_id <- bind_rows(hs %>% dplyr::rename(id = hs_gene_id),
							mm %>% dplyr::rename(id = mm_gene_id),
							mf %>% dplyr::rename(id = mf_gene_id),
							gg %>% dplyr::rename(id = gg_gene_id))

dbWriteTable(scEiaD, 'gene_name_id', name_id, overwrite = TRUE)
#dbWriteTable(scEiaD, 'haystack', haystack, overwrite = TRUE)
#db_create_index(scEiaD, table = 'haystack', columns = c('Gene'))

system2('touch', args[4])
