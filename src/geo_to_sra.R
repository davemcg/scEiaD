# Rscript

## extract sra project accessions with GEO ID in the dumbest way possible ->
## by building a url and dumping the web page into R

library(httr)
library(tidyverse)
library(DBI)

# web search to get info on GEO IDs
base_url <- 'https://www.ncbi.nlm.nih.gov/gds/?term='
GEO_ids <- scan('data/GEO_IDs.txt', what = 'character')
search_url <- paste0(base_url, paste(GEO_ids, collapse = '+'))

gse_prj <- data.frame(matrix(ncol = 2, nrow = 0))
# extract GEO ID and SRA project ID from web page
colnames(gse_prj) <- c("GSE", "SRA_PROJECT_ID")
for (i in GEO_ids){
  search_url <- paste0(base_url, i)
  page <- GET(search_url) %>% content(., 'text', )
  split <- (str_split(page, pattern = ' '))[[1]]
  acc <- split[grepl('acc', split)]
  acc <- (str_split(acc, pattern = '\\=|\\"|\\[|\\]|\\<|\\>')) %>% unlist()
  gse <- acc[grepl('GSE', acc)] %>% unique() 
  prj <- acc[grepl('PRJNA', acc)] %>% unique() 
  if (length(prj) == 0){prj <- NA}
  line = c(gse, prj)
  names(line) <- c('GSE', 'SRA_PROJECT_ID')
  gse_prj <- bind_rows(gse_prj, line)
}
# if SRA ID is missing, then load specific GEO page from GSE
# and search for it
base_url <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
missing_sra_gse <- gse_prj %>% filter(is.na(SRA_PROJECT_ID)) %>% pull(GSE)
gse_prj2 <- data.frame(matrix(ncol = 2, nrow = 0))
# extract GEO ID and SRA project ID from web page
colnames(gse_prj2) <- c("GSE", "SRA_PROJECT_ID")
for (i in missing_sra_gse){
  search_url <- paste0(base_url, i)
  page <- GET(search_url) %>% content(., 'text', )
  split <- (str_split(page, pattern = ' '))[[1]]
  acc <- split[grepl('acc', split)]
  acc <- (str_split(acc, pattern = '\\=|\\"|\\[|\\]|\\<|\\>')) %>% unlist()
  gse <- i
  prj <- acc[grepl('PRJNA', acc)] %>% unique() 
  if (length(prj) == 0){prj <- NA}
  line = c(gse, paste(prj, collapse = ','))
  names(line) <- c('GSE', 'SRA_PROJECT_ID')
  gse_prj2 <- bind_rows(gse_prj2, line)
}

library(GEOquery)
x <- getGEO('GSE12601')


# do not select all columns, as that costs more money
# bigquery essentially charges by data searched
# and if you exclude columns you reduce the amount
# of data you search across
con %>% 
  tbl('sra_study') %>% 
  select(accession, title, study_accession) %>% 
  filter(accession == 'SRP006387')