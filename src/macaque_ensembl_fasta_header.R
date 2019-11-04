args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)

headers <-read_lines(args[1])
ensmfat <- str_split(x, ' ') %>% map(1) %>% unlist()
gene <- str_split(x, ' ') %>% map(4) %>% unlist() %>% gsub('gene:', '', .)
symbol <- x %>% str_extract(., 'gene_symbol:\\S*') %>% gsub('gene_symbol:','',.)

out <- cbind(ensmfat, gene, symbol) %>% as_tibble()
write_tsv(out, path = args[2], col_names = FALSE)