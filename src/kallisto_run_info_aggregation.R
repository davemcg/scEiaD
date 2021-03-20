library(tidyverse)

dirs <- list.files('/data/OGVFB_BG/scEiaD/2021_01_21/', recursive=TRUE, pattern = 'run_info.json', full.names = TRUE)

run_info <- list()

for (i in dirs){run_info[[i]] <- jsonlite::fromJSON(file = i) %>% as_tibble()}

run_info <- run_info %>% bind_rows(.id = 'file') %>% mutate(sample_accession = str_extract(file, '(SRX|EGAF|ERS|SRS|iPSC_RPE_scRNA_)\\d+'))

write_tsv(run_info, file  = 'aggregated_run_info.tsv.gz')
