library(tidyverse)

f <- list.files('quant/', recursive = TRUE, pattern = 'run_info.json', full.names=TRUE)

f2 <- grep('old', f, value = TRUE, invert = TRUE)

run_info <- list()

for (i in f2){run_info[[i]] <- jsonlite::fromJSON(i) %>% as_tibble()}

run_info <- run_info %>% bind_rows(.id = 'file') %>% mutate(sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+'))

write_tsv(run_info, path = 'quant/aggregated_run_info.tsv.gz')
