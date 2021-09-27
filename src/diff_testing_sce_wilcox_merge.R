library(tidyverse)
library(scran)


make_data_tidy <- function(scran_diff_obj, group){
	load(scran_diff_obj)
	print(scran_diff_obj)
	tidy_data <- list()
	for (i in names(markers_wilcox)){
		if (group == 'cluster'){
			group <- 'Cluster'
		} 
		if (group == 'CellType_predict'){
			group <- 'CellType (Predict)'
		} 
		print(group)	
   		tidy_data[[i]] <- markers_wilcox[i] %>% as_tibble(rownames = 'Gene') %>%
   	    	mutate(Base = i) %>%
   	     	pivot_longer(cols = starts_with('AUC'),  names_to = 'Tested Against', values_to = 'AUC')
	}
	bind_rows(tidy_data) %>%
                    dplyr::select(-group, -group_name) %>%
                    mutate(Group = group,
							`Tested Against` = gsub('\\.',' ', `Tested Against`) %>% gsub('AUC ','',.),
                            `Tested Against` = case_when(`Tested Against` == 'AC HC_Precurs' ~ 'AC/HC Precursors',
                                                            TRUE ~ `Tested Against`))

}

args = commandArgs(trailingOnly=TRUE)
output <- args[length(args)]
group <- args[length(args) - 1]
input <- args[-length(args)]
input <- input[-length(input)]

tidy_list <- list()
for (i in input){
	tidy_list[[i]] <- make_data_tidy(i, group)
}

diff_testing <- tidy_list %>% bind_rows()
save(diff_testing, file = output)
