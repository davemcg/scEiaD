library(tidyverse)
library(glue)
library(rvest)

# get all SRP\d+
srp <- meta %>% 
  filter(grepl('10x', Platform), grepl('^SR|ER', sample_accession), !grepl('MTAB', study_accession)) %>% 
  pull(study_accession) %>% 
  gsub('#','',.) %>% 
  unique()

# harvest ncbi web page and get bam downloads
# if they don't exist, then these are files we got from SRA-dump in fastq format
# probably
bam_info<- list()
for (i in srp){
  print(i)
  ncbi_html <- read_html(glue("https://trace.ncbi.nlm.nih.gov/Traces/sra/?study={i}")) 
  bam_name <- ncbi_html %>% 
    html_nodes('a') %>% 
    html_text() %>% 
    grep('bam',., value = TRUE)
  bam_url <- ncbi_html %>% 
    html_nodes('a') %>% 
    grep('bam',., value = TRUE) %>% 
    str_extract(., 'http.*\\"') %>% 
    gsub('\\"','',.)
  out <- bind_cols(bam_name,bam_url)
  colnames(out) <- c('bam','url')
  bam_info[[i]] <- out
  Sys.sleep(1)
}

# now get SRP <-> SRS <-> SRX <-> SRR mapping
# "hand" search ncbi with `paste(srp, collapse = ' OR ')`
# save as -> file -> xml
library(xml2)
full_xml <- read_xml('~/Downloads/SraExperimentPackage.xml') %>% as_list()
exp_set <- full_xml$EXPERIMENT_PACKAGE_SET
exp_set_list <- lapply(exp_set, function(x) unlist(x, recursive = T))
common_names <- lapply(exp_set_list, names) %>% reduce(intersect)
res <-lapply(exp_set_list , function(x) x[common_names]) %>%   do.call(rbind,.) %>% as.data.frame()

sra_meta <- res %>% select('EXPERIMENT.IDENTIFIERS.PRIMARY_ID', 'RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID', 'Pool.Member.IDENTIFIERS.PRIMARY_ID', 'EXPERIMENT.STUDY_REF.IDENTIFIERS.PRIMARY_ID') %>% as_tibble()
colnames(sra_meta) <- c('sample_accession_2', 'run_accession', 'sample_accession', 'study_accession')



our_meta <- read_tsv('~/git/scEiaD/data/sample_run_layout_organism_tech.tsv')

a <- our_meta %>% 
  select(sample_accession) %>% 
  unique() %>% 
  left_join(
    bind_rows(bam_info, .id = 'study_accession') %>% 
      mutate(run_accession = str_extract(url, 'SRR\\d+')) %>% 
      left_join(sra_meta) %>% 
      select(-sample_accession_2),
    by = 'sample_accession') %>% 
  filter(!is.na(bam))

b <- our_meta %>% 
  select(sample_accession) %>% 
  unique() %>% 
  left_join(
    bind_rows(bam_info, .id = 'study_accession') %>% 
                mutate(run_accession = str_extract(url, 'SRR\\d+')) %>% 
                left_join(sra_meta) %>% 
      select(-sample_accession),
              by = c('sample_accession' = 'sample_accession_2')) %>% 
  filter(!is.na(bam))

download_table <- bind_rows(a,b)
write_tsv(download_table, file = 'data/sra_10x_bam_downloads.tsv')
