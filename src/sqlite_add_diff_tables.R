library(tidyverse)
library(pool)
library(RSQLite)

args = commandArgs(trailingOnly=TRUE)


anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)

prefix1 = args[2]
prefix2 = args[3]

# A
# 'Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__200__0.1__7__100.A.diff.coef_table.Rdata'
# 'Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__dims200__7__100__0.1__cluster.sceWilcox_summary.Rdata'
load(paste0(prefix1, '.A.diff.coef_table.Rdata'))
load(paste0(prefix2, '__cluster.sceWilcox_summary.Rdata'))
coefficient_table <- coefficient_table %>% filter(grepl('seuratCluster', term)) %>% mutate(term = gsub('seuratCluster','',term))
coefficient_table <- coefficient_table %>% left_join(markers_summary, by = c('term' = 'cluster', 'gene_short_name' = 'gene'))
dbWriteTable(anthology_2020_v01, 'diff_testing_A', coefficient_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'diff_testing_A', columns = c('gene_short_name'))
db_create_index(anthology_2020_v01, table = 'diff_testing_A', columns = c('term'))


# C
# 'Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__200__0.1__7__100.C.diff.coef_table.Rdata'
# 'Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__dims200__7__100__0.1__CellType_predict.sceWilcox_summary.Rdata'
load(paste0(prefix1, '.C.diff.coef_table.Rdata'))
load(paste0(prefix2, '__CellType_predict.sceWilcox_summary.Rdata'))
coefficient_table <- coefficient_table %>% filter(grepl('CellType_predict', term)) %>% mutate(term = gsub('CellType_predict','',term))
coefficient_table <- coefficient_table %>% left_join(markers_summary, by = c('term' = 'cluster', 'gene_short_name' = 'gene'))
dbWriteTable(anthology_2020_v01, 'diff_testing_C', coefficient_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'diff_testing_C', columns = c('gene_short_name'))
db_create_index(anthology_2020_v01, table = 'diff_testing_C', columns = c('term'))

# E
load(paste0(prefix1, '.E.diff.coef_table.Rdata'))
load(paste0(prefix2, '__CellType.sceWilcox_summary.Rdata'))
coefficient_table <- coefficient_table %>% filter(grepl('CellType', term)) %>% mutate(term = gsub('CellType','',term))
coefficient_table <- coefficient_table %>% left_join(markers_summary, by = c('term' = 'cluster', 'gene_short_name' = 'gene'))
dbWriteTable(anthology_2020_v01, 'diff_testing_E', coefficient_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'diff_testing_E', columns = c('gene_short_name'))
db_create_index(anthology_2020_v01, table = 'diff_testing_E', columns = c('term'))

# G 
load(paste0(prefix1, '.G.SC.diff.coef_table.Rdata'))
load(paste0(prefix2, '__subcluster.sceWilcox_summary.Rdata'))
coefficient_table <- coefficient_table %>% filter(grepl('seuratSubCluster', term)) %>% mutate(term = gsub('seuratSubCluster','',term))
coefficient_table <- coefficient_table %>% left_join(markers_summary, by = c('term' = 'cluster', 'gene_short_name' = 'gene'))
dbWriteTable(anthology_2020_v01, 'diff_testing_G', coefficient_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'diff_testing_G', columns = c('gene_short_name'))
db_create_index(anthology_2020_v01, table = 'diff_testing_G', columns = c('term'))

# doublets
load(args[4])
dbWriteTable(anthology_2020_v01, 'doublets', doublet_call_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'doublets', columns = c('Barcode'))
