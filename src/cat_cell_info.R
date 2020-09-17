library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input <- list()
for (i in args){
	input[[i]] <- read_tsv(i)
}

out <- input %>% bind_rows()

write_tsv(out, path = 'cell_info.tsv')
