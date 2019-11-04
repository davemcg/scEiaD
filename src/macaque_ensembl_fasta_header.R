args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)

headers <-read_lines(args[1])
ensmfat <- str_split(headers, ' ') %>% map(1) %>% unlist()
gene <- str_split(headers, ' ') %>% map(4) %>% unlist() %>% gsub('gene:', '', .)
symbol <- headers %>% str_extract(., 'gene_symbol:\\S*') %>% gsub('gene_symbol:','',.)

out <- cbind(ensmfat, gene, symbol) %>% as_tibble() %>% 
  mutate(symbol = case_when(is.na(symbol) ~ gene,
                            TRUE ~ symbol))
write_tsv(out, path = args[2], col_names = FALSE)