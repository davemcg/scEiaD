library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
meta <- read_tsv(args[1])

runs <- meta$run_accession
exists <- list.files(args[2])

missing <- list()
for (i in runs){
	hits = grep(i, exists)
	if (length(hits) == 0){
		missing[[i]] <- i
	}
}

write_tsv(missing %>% unlist() %>% as_tibble(), file = args[3])
