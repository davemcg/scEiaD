library(tidyverse)
library(xml2)
library(rvest)
args = commandArgs(trailingOnly=TRUE)
# example args:
# ~/Downloads/SraExperimentPackage.xml (# https://www.ncbi.nlm.nih.gov/sra/?term=SRP238072 -> Send to -> File -> Full XML)
# SRP238072
# SRP238072_meta.tsv.gz

# full_xml <- read_xml('~/Downloads/SraExperimentPackage.xml') %>% as_list()
full_xml <- read_xml(args[1]) %>% as_list()
exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
common_names <- lapply(exp_set_list, names) %>% reduce(intersect)
res <-lapply(exp_set_list , function(x) x[common_names]) %>% do.call(rbind,.) %>% as.data.frame()
# write_tsv(res, '~/Downloads/SRP286543.tsv')

# get bam url
url <-  paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=", args[2])
links <- url %>% 
  read_html() %>% 
  html_nodes("td:nth-child(3) a") %>% 
  html_attr("href") %>% 
  as_tibble() %>% 
  mutate(`RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID` = str_extract(value, 'SRR\\d+' )) %>% 
  rename(bam_url = value)

res <- res %>% left_join(links)
write_tsv(res, file = args[3])