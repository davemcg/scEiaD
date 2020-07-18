library(tidyverse)
library(pool)
library(RSQLite)

args = commandArgs(trailingOnly=TRUE)


anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)

load(args[2])
load(args[3])

PB_results <- bind_rows(PB_resultsC2, PB_resultsABC)

dbWriteTable(anthology_2020_v01, 'PB_results', PB_results, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'PB_results', columns = c('Gene'))
db_create_index(anthology_2020_v01, table = 'PB_results', columns = c('test'))
db_create_index(anthology_2020_v01, table = 'PB_results', columns = c('PB_Test'))

# doublets
load(args[4])
dbWriteTable(anthology_2020_v01, 'doublets', doublet_call_table, overwrite = TRUE)
db_create_index(anthology_2020_v01, table = 'doublets', columns = c('Barcode'))
