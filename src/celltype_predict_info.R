args = commandArgs(trailingOnly=TRUE)

library(tidyverse)



grab_var <- function(feature, obj_path){
	str_extract(obj_path, paste0(feature, '-[\\d|\\w]+_')) %>% gsub(paste0(feature,'-|__'),'',.)
}

calculator <- function(obj_path){
	load(obj_path)
	vars <- list()
	vars$n_features <- grab_var('n_features', obj_path)
	vars$transform <- grab_var('transform', obj_path)
	vars$partition <- grab_var('partition', obj_path)
	vars$covariate <- grab_var('covariate', obj_path)
	vars$method <- grab_var('method', obj_path)
	vars$dims <- grab_var('dims', obj_path)
	vars$epochs <- grab_var('epochs', obj_path)
	# celltype_predict counts
	out <- list()
	out$CTP_counts <- umap$CellType_predict %>% table() %>% enframe(name = 'CellType_predict', value = 'Count')
	out$CTP_counts <- bind_cols(out$CTP_counts, vars %>% unlist() %>% rbind() %>% as_tibble())
	out$CTP_organism_counts <- umap %>% group_by(CellType_predict, organism) %>% summarise(Count = n())
	out$CTP_organism_counts <- bind_cols(out$CTP_organism_counts, vars %>% unlist() %>% rbind() %>% as_tibble())
	out$CTP_study_counts <- umap %>% group_by(CellType_predict, study_accession) %>% summarise(Count = n())
	out$CTP_study_counts <- bind_cols(out$CTP_study_counts, vars %>% unlist() %>% rbind() %>% as_tibble())
	out$bharti_CTP_counts <-  umap %>% filter(grepl('Bharti', study_accession)) %>% group_by(CellType_predict) %>% summarise(Count = n()) %>% arrange(-Count)
	out$bharti_CTP_counts <- bind_cols(out$bharti_CTP_counts, vars %>% unlist() %>% rbind() %>% as_tibble())
	out$na_count <- umap %>% filter(is.na(CellType_predict) & !is.na(TabulaMurisCellType_predict)) %>% nrow()
	out
}

data <- list()
files <- args[-length(args)]
output <- args[length(args)]
for (i in files){
	data[[i]] <- calculator(i)
}

# bind tibbles together
CTP_counts <- data %>% map('CTP_counts') %>% bind_rows()
CTP_organism_counts <- data %>% map("CTP_organism_counts") %>% bind_rows()
CTP_study_counts <- data %>% map("CTP_study_counts") %>% bind_rows()
bharti_CTP_counts <- data %>% map("bharti_CTP_counts") %>% bind_rows()

CTP_data_list <- data
save(CTP_counts, CTP_organism_counts, CTP_study_counts, bharti_CTP_counts, CTP_data_list, file = output)
# bharti_CTP_counts %>% group_by(n_features, transform, partition, covariate, method, dims) %>% mutate(Ratio = Count / sum(Count)) %>% filter(CellType_predict == 'RPE') %>% arrange(-Ratio) %>% ungroup() %>% select( CellType_predict, Count, n_features, dims, epochs, Ratio) %>% data.frame()
# bharti_CTP_counts %>% group_by(n_features, transform, partition, covariate, method, dims) %>% mutate(Ratio = Count / sum(Count)) %>% filter(CellType_predict == 'Muller Glia') %>% arrange(-Ratio) %>% ungroup() %>% select( CellType_predict, Count, n_features, dims, epochs, Ratio) %>% data.frame()
