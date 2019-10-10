library(tidyverse)

meta <- read_tsv('data/sample_run_layout_organism_tech.tsv')
fastq_path <- '/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/fastq/'
fastq_files_existing <- list.files(fastq_path, pattern = '*.fastq.gz')

fastq_wanted <- meta %>% 
  filter(grepl('^SRS', sample_accession)) %>% 
  mutate(fastqS = case_when(library_layout == 'SINGLE' ~ paste0(run_accession, '.fastq.gz')),
         fastq1 = case_when(library_layout == 'PAIRED' ~ paste0(run_accession, '_1.fastq.gz')),
         fastq2 = case_when(library_layout == 'PAIRED' ~ paste0(run_accession, '_2.fastq.gz'))) %>% 
  select(run_accession, fastqS, fastq1, fastq2) %>% 
  group_by(run_accession) %>% 
  summarise(fastq = list(c(fastqS, fastq1, fastq1))) %>% 
  unnest(fastq) %>% 
  pull(fastq) %>% 
  unique()

fastq_wanted <- fastq_wanted[!is.na(fastq_wanted)]

(fastq_wanted %in% fastq_files_existing) %>% table()

fastq_missing <- fastq_wanted[!(fastq_wanted %in% fastq_files_existing)]
srr <- fastq_missing %>% gsub('_2.fastq.gz|_1.fastq.gz|.fastq.gz', '', .) %>% unique()

swarm <- srr %>% 
  enframe() %>% 
  glue::glue_data('fasterq-dump --include-technical --split-files {value}; pigz -f -p 8 {value}*fastq; rm {value}*fastq')
write(swarm, file = 'data/fasterq_dump_2019_10_10.swarm')
