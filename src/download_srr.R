library(tidyverse)
load('data/sra_metadata.Rdata')
# write swarm fasterq dump command for biowulf
# cd /data/mcgaugheyd/projects/nei/mcgaughey/massive_integrated_eye_scRNA/fastq
# swarm -f ~/git/massive_integrated_eye_scRNA/data/fasterq_dump.swarm -b 400 --time=0:30:00 --module=sratoolkit -t 4 -g 16
write(sra_metadata %>% mutate(fastq_command = paste0('fasterq-dump ', run_accession, '; pigz -p 4 ', run_accession, '*fastq')) %>% pull(fastq_command), file = 'src/fasterq_dump.swarm')
# existing swarm command which removes already downlaoded files
existing <- list.files(path = '/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/fastq/', pattern = "fastq.gz") %>% gsub('.fastq.gz','',.) %>% gsub('_.*','',.)
write(sra_metadata %>% filter(!run_accession %in% existing) %>% mutate(fastq_command = paste0('fasterq-dump ', run_accession, '; pigz -p 4 ', run_accession, '*fastq')) %>% pull(fastq_command), file = 'src/fasterq_dump_missing.swarm')


# ID paired / single end discrepancy
# paired end supposed to have _1 and _2
existing <- list.files(path = '/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/fastq/', pattern = "fastq.gz") 

files = ''
for (idx in seq(1,nrow(sra_metadata))){
  if (sra_metadata[idx,'library_layout'] == 'PAIRED'){
    files = c(files, paste0(sra_metadata[idx,'run_accession'], '_1.fastq.gz'),
              paste0(sra_metadata[idx,'run_accession'], '_2.fastq.gz'))
  } else {
    files = c(files, paste0(sra_metadata[idx,'run_accession'], '.fastq.gz'))
  }
}
files <- files[-1]
sra_metadata %>% filter(run_accession %in% (files[(!files %in% existing)] %>% gsub('_.*','',.))) %>% 
  mutate(fastq_command = paste0('fasterq-dump --include-technical --split-files', run_accession, '; pigz -p 8 ', run_accession, '*fastq')) %>% pull(fastq_command)
