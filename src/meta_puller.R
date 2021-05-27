library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(xml2)
library(rvest)
args = commandArgs(trailingOnly=TRUE)

# example args:
# SRP238072
# SRP238072_meta.tsv.gz

xml_path <- paste0('http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=fullxml&term=', args[1], '[Experiment]')
# grab sra xml data for accession
full_xml <- read_xml(xml_path) %>% as_list()
exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
common_names <- lapply(exp_set_list, names) %>% reduce(intersect)
res <-lapply(exp_set_list , function(x) x[common_names]) %>% do.call(rbind,.) %>% as.data.frame()

# get sra "RunInfo"
runinfo_path <- paste0('http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=', args[1], '[Experiment]')
runinfo <- read_csv(runinfo_path)
res <- res %>% left_join(., runinfo, by = c("RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID" = "Run"))

# get bam url
url <-  paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=", args[1])
links <- url %>% 
  read_html() %>% 
  html_nodes("td:nth-child(3) a") %>% 
  html_attr("href") %>% 
  as_tibble() %>% 
  mutate(`RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID` = str_extract(value, 'SRR\\d+' )) %>% 
  rename(bam_url = value)

res <- res %>% left_join(links, by = "RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID")
write_tsv(res, file = args[2])