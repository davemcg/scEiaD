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
diff_testing_CTp <- diff_testing_wilcox
diff_testing_CTp_summary <- diff_summary_wilcox

load(args[3]) # CellType
diff_testing_CT <- diff_testing_wilcox
diff_testing_CT_summary <- diff_summary_wilcox

load(args[4]) # cluster
diff_testing_cluster <- diff_testing_wilcox
diff_testing_cluster_summary <- diff_summary_wilcox

wilcox_diff_testing <- bind_rows(diff_testing_CTp, diff_testing_CT, diff_testing_cluster) %>% arrange(-AUC)
wilcox_diff_testing_geneVals <- bind_rows(diff_testing_CTp_summary, diff_testing_CT_summary, diff_testing_cluster_summary) %>% dplyr::rename(Gene = gene, p.value = pval)

# yank out pval and FDR info (which are at the base - gene level)
wilcox_diff_AUC <- wilcox_diff_testing %>% select(Gene, Base, `Tested Against`, AUC, logFC, Group)


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
wilcox_diff_AUC_human <- wilcox_diff_AUC %>% filter(grepl('ENSG', Gene)) %>% left_join(gene_id_converter %>% dplyr::select(hs_gene_id, hs_gene_name, human_name) %>% unique(), by = c('Gene' = 'hs_gene_id')) %>% mutate(Gene = paste0(hs_gene_name, ' (', Gene, ')'))  %>% dplyr::select(-hs_gene_name)

wilcox_diff_AUC_mouse <- wilcox_diff_AUC %>% filter(grepl('ENSMUS', Gene)) %>% left_join(gene_id_converter %>% dplyr::select(mm_gene_id, mm_gene_name, human_name) %>% unique(), by = c('Gene' = 'mm_gene_id')) %>% mutate(Gene = paste0(mm_gene_name, ' (', Gene, ')'))  %>% dplyr::select(-mm_gene_name)

wilcox_diff_AUC <- bind_rows(wilcox_diff_AUC_human, wilcox_diff_AUC_mouse) %>% dplyr::rename(`Gene Name` = human_name) %>% mutate(AUC = as.numeric(AUC)) %>% arrange(-AUC)
dbWriteTable(scEiaD, 'wilcox_diff_AUC', wilcox_diff_AUC, overwrite = TRUE)
db_create_index(scEiaD, table = 'wilcox_diff_AUC', columns = c('Gene'))
db_create_index(scEiaD, table = 'wilcox_diff_AUC', columns = c('Base'))
db_create_index(scEiaD, table = 'wilcox_diff_AUC', columns = c('Group'))
db_create_index(scEiaD, table = 'wilcox_diff_AUC', columns = c('Tested Against'))
wilcox_group_base_sets <- wilcox_diff_AUC %>% dplyr::select(Group, Base)  %>% unique() %>% arrange(Group, Base)
dbWriteTable(scEiaD, 'wilcox_diff_AUC_sets', wilcox_group_base_sets, overwrite = TRUE)
wilcox_diff_AUC_genes <- wilcox_diff_AUC %>% pull(Gene) %>% unique()
wilcox_diff_AUC_genes <- wilcox_diff_AUC_genes %>% tibble::enframe(value = 'Gene') %>% dplyr::select(-name) %>% arrange()
dbWriteTable(scEiaD, 'wilcox_diff_AUC_genes', wilcox_diff_AUC_genes, overwrite = TRUE)
## now the pval table
wilcox_diff_testing_geneVals__human <- wilcox_diff_testing_geneVals %>% filter(grepl('ENSG', Gene)) %>% left_join(gene_id_converter %>% dplyr::select(hs_gene_id, hs_gene_name, human_name) %>% unique(), by = c('Gene' = 'hs_gene_id')) %>% mutate(Gene = paste0(hs_gene_name, ' (', Gene, ')'))  %>% dplyr::select(-hs_gene_name)
wilcox_diff_testing_geneVals__mouse <- wilcox_diff_testing_geneVals  %>% filter(grepl('ENSMUS', Gene)) %>% left_join(gene_id_converter %>% dplyr::select(mm_gene_id, mm_gene_name, human_name) %>% unique(), by = c('Gene' = 'mm_gene_id')) %>% mutate(Gene = paste0(mm_gene_name, ' (', Gene, ')'))  %>% dplyr::select(-mm_gene_name)
wilcox_diff_testing <- bind_rows(wilcox_diff_testing_geneVals__human, wilcox_diff_testing_geneVals__mouse) %>% dplyr::rename(`Gene Name` = human_name) %>% 
							mutate(p.value = as.numeric(p.value), mean_auc = as.numeric(mean_auc)) %>%
							arrange(`p.value`, -mean_auc) %>% dplyr::rename(Base = cluster, auc_count_070 = count)
dbWriteTable(scEiaD, 'wilcox_diff_testing', wilcox_diff_testing, overwrite = TRUE)
db_create_index(scEiaD, table = 'wilcox_diff_testing', columns = c('Gene'))
db_create_index(scEiaD, table = 'wilcox_diff_testing', columns = c('Base'))
db_create_index(scEiaD, table = 'wilcox_diff_testing', columns = c('Group'))




# doublets
print('write doublets')
load(args[5])
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

system2('touch', args[6])
