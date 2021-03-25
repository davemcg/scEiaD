library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

accuracy_list <- list()
for (i in args){
	load(i)
	accuracy_list[[i]] <- accuracy %>% mutate(file = i)
}

accuracy <- bind_rows(accuracy_list) 
save(accuracy, file = 'pipeline_data/results/merged_xgboost.Rdata')
 
