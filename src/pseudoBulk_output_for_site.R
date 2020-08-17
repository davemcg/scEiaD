args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(edgeR)

load(args[1])
####### pseudobulk Objs############33

# fix column naming
study = gsub('donor_','donor', edgeR_obj$processed_data$cts %>% colnames) %>% str_split(., '_') %>% map(., 1)
platform = gsub('donor_','donor', edgeR_obj$processed_data$cts %>% colnames) %>% str_split(., '_') %>% map(., 2)
covariate = gsub('donor_','donor', edgeR_obj$processed_data$cts %>% colnames) %>% str_split(., '_') %>% map(., 3)
term = gsub('donor_','donor', edgeR_obj$processed_data$cts %>% colnames) %>% str_split(., '_') %>% map(., 4)
organism = gsub('donor_','donor', edgeR_obj$processed_data$cts %>% colnames) %>% str_split(., '_') %>% map(., 5)

pseudoBulk = edgeR_obj$processed_data$cts %>% as_tibble()
colnames(pseudoBulk) = paste(study, platform, covariate, term, organism, sep = '__')

write_tsv(pseudoBulk, path = args[2])
