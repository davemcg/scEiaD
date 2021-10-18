library(tidyverse)
library(scran)


make_data_tidy <- function(scran_diff_obj, group, test = 'wilcox'){
	load(scran_diff_obj)
	print(scran_diff_obj)
	tidy_data <- list()
	if (test == 'wilcox'){
		col_val <- 'AUC'
	} else {
		col_val <- 'logFC'
	}
	for (i in names(markers_wilcox)){
		if (group == 'cluster'){
			group <- 'Cluster'
		} 
		if (group == 'CellType_predict'){
			group <- 'CellType (Predict)'
		}
		if (test == 'wilcox'){
   			tidy_data[[i]] <- markers_wilcox[i] %>% as_tibble(rownames = 'Gene') %>%
   	    		mutate(Base = i) %>%
   	     		pivot_longer(cols = starts_with(col_val),  names_to = 'Tested Against', values_to = col_val)
		} else if (test == 't'){
   			tidy_data[[i]] <- markers_t[i] %>% as_tibble(rownames = 'Gene') %>%
   	    		mutate(Base = i) %>%
   	     		pivot_longer(cols = starts_with(col_val),  names_to = 'Tested Against', values_to = col_val)
		}
	}
	bind_rows(tidy_data) %>%
                    dplyr::select(-group, -group_name) %>%
                    mutate(Group = group,
							`Tested Against` = gsub('\\.',' ', `Tested Against`) %>% gsub(paste(col_val, ''),'',.),
                            `Tested Against` = case_when(`Tested Against` == 'AC HC_Precurs' ~ 'AC/HC Precursors',
                                                            TRUE ~ `Tested Against`))

}

args = commandArgs(trailingOnly=TRUE)
output <- args[length(args)]
group <- args[length(args) - 1]
input <- args[-length(args)]
input <- input[-length(input)]

tidy_list_wilcox <- list()
for (i in input){
	tidy_list_wilcox[[i]] <- make_data_tidy(i, group, 'wilcox')
}
tidy_list_t <- list()
for (i in input){
	tidy_list_t[[i]] <- make_data_tidy(i, group, 't')
}

# process summary data
input_summary <- gsub('.Rdata', '_summary.Rdata', input) 
summary_list_wilcox <- list()
for (i in input_summary){
	print(i)
	load(i)
	markers_summary_w$Group = group
    summary_list_wilcox[[i]] <- markers_summary_w
}
summary_list_t <- list()
for (i in input_summary){
	print(i)
	load(i)
	markers_summary_t$Group = group
    summary_list_t[[i]] <- markers_summary_t
}



diff_summary_wilcox <- summary_list_wilcox %>% bind_rows() %>% mutate(Group = case_when(Group == 'CellType_predict' ~ 'CellType (Predict)',
																				Group == 'cluster' ~ 'Cluster',
																				TRUE ~ Group))
diff_testing_wilcox <- tidy_list_wilcox %>% bind_rows()

diff_summary_t <- summary_list_t %>% bind_rows() %>% mutate(Group = case_when(Group == 'CellType_predict' ~ 'CellType (Predict)',
                                                                                Group == 'cluster' ~ 'Cluster',
                                                                                TRUE ~ Group))
diff_testing_t <- tidy_list_t %>% bind_rows()


# add logFC to wilcox data 
diff_summary_wilcox <- diff_summary_wilcox %>% left_join(diff_summary_t %>% select(gene, mean_logFC, cluster, Group))
diff_testing_wilcox <- diff_testing_wilcox %>% left_join(diff_testing_t %>% select(Gene, Base, `Tested Against`, Group, logFC))

save(diff_testing_wilcox, diff_summary_wilcox, diff_testing_t, diff_summary_t, file = output)
