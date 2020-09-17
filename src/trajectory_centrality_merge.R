library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

input <- args[-1]
out <- args[length(args)]

data <- list()
for (i in input){
	data[[i]] <- read_tsv(i)
}

save(data, file = out)
