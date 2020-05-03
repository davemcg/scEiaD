library(tidyverse)
library(pool)
library(RSQLite)


anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-2000-counts-onlyDROPLET-batch-scVI-200-0.1-100-7.sqlite", idleTimeout = 3600000)



# A
load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__200__0.1__7__100.A.diff.coef_table.Rdata')
load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__dims200__7__100__0.1__cluster.sceWilcox_summary.Rdata')
coefficient_table <- coefficient_table %>% filter(grepl('seuratCluster', term)) %>% mutate(term = gsub('seuratCluster','',term))
coefficient_table <- coefficient_table %>% left_join(markers_summary, by = c('term' = 'cluster', 'gene_short_name' = 'gene'))
dbWriteTable(anthology_2020_v01, 'diff_testing_A', coefficient_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'diff_testing_A', columns = c('gene_short_name'))
db_create_index(anthology_2020_v01, table = 'diff_testing_A', columns = c('term'))


# C
load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__200__0.1__7__100.C.diff.coef_table.Rdata')
load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__2000__counts__onlyDROPLET__batch__scVI__dims200__7__100__0.1__CellType_predict.sceWilcox_summary.Rdata')
coefficient_table <- coefficient_table %>% filter(grepl('CellType_predict', term)) %>% mutate(term = gsub('CellType_predict','',term))
coefficient_table <- coefficient_table %>% left_join(markers_summary, by = c('term' = 'cluster', 'gene_short_name' = 'gene'))
dbWriteTable(anthology_2020_v01, 'diff_testing_C', coefficient_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'diff_testing_C', columns = c('gene_short_name'))
db_create_index(anthology_2020_v01, table = 'diff_testing_C', columns = c('term'))

