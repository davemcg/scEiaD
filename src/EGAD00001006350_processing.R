# scan EGAD00001006350 and extract metadata from this spaghetti structure
library(tidyverse)
library(zellkonverter)

f <- list.files('~/Downloads/EGAD00001006350/xmls/samples/', full.names = TRUE)

sample_meta <- list()
for (i in f){
  raw <- scan(i, what = 'character') %>% paste(., collapse = ' ')
  samp_name <- raw %>% str_extract(., '<TITLE>.*</TITLE>') %>% gsub('<TITLE>|</TITLE>', '', .)
  samp_info <- raw %>% str_extract(., '<DESCRIPTION>.*</DESCRIPTION>')%>% gsub('<DESCRIPTION>|</DESCRIPTION>', '', .)
  sample_meta[[i]] <- data.frame(sample = samp_name, info = samp_info) %>% as_tibble()
}
sample_meta <- sample_meta %>% bind_rows()
sample_meta %>% mutate()


o <- readH5AD('~/Downloads/adata_final_clean_organoid.h5ad')
f <- readH5AD('~/Downloads/adata_final_clean_fovea.h5ad')
p <- readH5AD('~/Downloads/adata_final_clean_periphery.h5ad')

coldata <- bind_rows(f@colData %>% as_tibble(rownames = 'row'), p@colData %>% as_tibble(rownames = 'row'), o@colData %>% as_tibble(rownames = 'row'))

coldata <- coldata %>% left_join(sample_meta, by = c('condition' = 'sample'))

map <- read_tsv('~/Downloads/EGAD00001006350/delimited_maps/Sample_File.map', col_names = FALSE)
colnames(map) <- c('ID', 'EGAN','bam','EGAF')

coldata <- coldata %>% mutate(ID = str_extract(info, '^.*_lib\\d+')) %>% left_join(map)

                      