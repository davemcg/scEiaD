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
colnames(gse_prj) <- c("GSE", "SRA_PROJECT_ID")
for (i in GEO_ids){
  search_url <- paste0(base_url, i)
  page <- GET(search_url) %>% content(., 'text', )
  split <- (str_split(page, pattern = ' |\\=|\\"|\\[|\\]|\\<|\\>'))[[1]]
  gse <- split[grepl('GSE',split)] %>% unique() %>% head(1)
  prj <- split[grepl('PRJNA',split)] %>% unique() %>% head(1)
  if (length(prj) == 0){prj <- NA}
  line = c(gse, prj)
  names(line) <- c('GSE', 'SRA_PROJECT_ID')

  gse_prj <- bind_rows(gse_prj, line)
}
page <- GET(search_url) %>% content(., 'text', )
split <- str_split(page, pattern = ' |\\=|\\"|\\[|\\]|\\<|\\>')
  # get SRA accessions / info from @seandavi12 omicidx 
  # billing is the project name in the Google
  # bigQuery web UI
  con <- dbConnect(bigrquery::bigquery(), 
                   project = 'isb-cgc-01-0006', 
                   dataset = 'omicdx', billing = 
                     'rosy-solstice-244919')

# do not select all columns, as that costs more money
# bigquery essentially charges by data searched
# and if you exclude columns you reduce the amount
# of data you search across
con %>% 
  tbl('sra_study') %>% 
  select(accession, title, study_accession) %>% 
  filter(accession == 'SRP006387')