library(tidyverse)
library(pool)
library(RSQLite)
library(Matrix)
args = commandArgs(trailingOnly=TRUE)
git_dir = Sys.getenv('SCIAD_GIT_DIR')
library(glue)
source(glue('{git_dir}/src/make_gene_id_converter_table.R'))

scEiaD <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)

load(args[2])
haystack <- scH$results %>% as_tibble(rownames = 'hs_gene_id') %>% left_join(gene_id_converter %>% dplyr::select(hs_gene_id, hs_gene_name) %>% unique()) %>% mutate(Gene = paste0(hs_gene_name, ' (', hs_gene_id, ')')) %>% dplyr::select(Gene, D_KL, `log.p.vals`, `log.p.adj`, `T.counts`)

dbWriteTable(scEiaD, 'haystack', haystack, overwrite = TRUE)
db_create_index(scEiaD, table = 'haystack', columns = c('Gene'))

# label most likely [ ] specificity for each gene
# where [ ] is clsuter, celtype, celltype_predict
ctp <- scEiaD %>% tbl('wilcox_diff_testing') %>% filter(Group == 'CellType (Predict)') %>% collect() %>% group_by(Gene) %>% slice_min(p.value,n = 3) %>% slice_max(mean_logFC, n = 1) %>% filter(FDR < 1 , mean_logFC > 0.1) %>% select(-`Gene Name`) %>% unique()


 allowed_ct <- scEiaD %>% tbl('wilcox_diff_testing') %>% filter(Group == 'CellType (Predict)') %>% collect() %>% pull(Base) %>% unique()
ct <- scEiaD %>% tbl('wilcox_diff_testing') %>% filter(Group == 'CellType', Base %in% allowed_ct) %>% collect() %>% group_by(Gene) %>% slice_min(p.value,n = 3) %>% slice_max(mean_logFC, n = 1) %>% filter(FDR < 1 , mean_logFC > 0.1) %>% select(-`Gene Name`) %>% unique()

cluster <- scEiaD %>% tbl('wilcox_diff_testing') %>% filter(Group == 'Cluster') %>% collect() %>% group_by(Gene) %>% slice_min(p.value,n = 3) %>% slice_max(mean_logFC, n = 1) %>% filter(FDR < 1, mean_logFC > 0.1) %>% select(-`Gene Name`) %>% unique()

gene_auto_label <- ctp %>% select(Gene, CellType_predict = Base) %>% full_join(cluster %>% select(Gene, Cluster = Base))  %>% full_join(ct %>% select(Gene, CellType = Base))

dbWriteTable(scEiaD, 'gene_auto_label', gene_auto_label, overwrite = TRUE)
db_create_index(scEiaD, table = 'gene_auto_label', columns = c('Gene'))
