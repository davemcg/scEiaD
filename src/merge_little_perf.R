library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

files <- list.files('sample_level/perf', full.names = TRUE)
files <- files[grepl('tsv', files)]
files <- files[grepl('norm', files)]

perf <- list()
for (i in files){
	perf[[i]] <- read_tsv(i)
}

perf <- bind_rows(perf)
write_tsv(perf, file = args[1])



