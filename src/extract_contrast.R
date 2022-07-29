args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(DESeq2)
library(parallel)

load(args[1])

group <- args[2] #  'CellType_predict'

group_vals <- colData(deseq2_obj)[,group] %>% as.character() %>% unique()
colData(deseq2_obj)$cluster <- as.factor(colData(deseq2_obj)$cluster)
get_contrast <- function(target, against = 'All'){
	
	# this function runs the target (e.g. Cone against all remaining data) 
	# with a contrast pull for DESeq2
	target_dot <- gsub('\\/|-| |\\(|\\)|\'', '.', target)
	# remove the target celltype/cluster
	group_vals_minus <- group_vals[!group_vals %in% target]
	# https://support.bioconductor.org/p/86347/
	if (against == 'All'){
		diff <- try({results(deseq2_obj,
   		           contrast = list(c(paste0(group, target_dot)), paste0(group, gsub('\\/|-| |\\(|\\)|\'', '.', group_vals_minus))),
   	    	       listValues = c(1,-1/(length(group_vals_minus) + 1 )))	
		})
		 diff <- try({diff %>% as_tibble(rownames = 'Gene') })
		 diff$Against <- against
	} else {
		# identify unique tests
		combos <- combn(group_vals, 2) %>% data.frame() %>% t()
		colnames(combos) <- c('a','b')
		combos <- combos %>% as_tibble() %>% mutate(pair = glue::glue("{a}_{b}"))
	
		diff <- list()
		for (i in group_vals_minus){
			combo <- paste0(target, '_', i)
			if (combo %in% combos$pair){
				print(paste(group, combo))

				diff[[i]] <- try({results(deseq2_obj,
   		           contrast = c(group, target, i))	
				})
			
				diff[[i]]$Against <- i
				diff[[i]] <- try({ diff[[i]]  %>% as_tibble(rownames = 'Gene') })
			}
		}	
		diff <- bind_rows(diff)
	}

	if (class(diff) != 'try-error'){
		diff$Group <- group
		diff$Base <- target
	}
	diff
}


de_results <- mclapply(group_vals, get_contrast, mc.cores = 16)
names(de_results) <- group_vals
de_results_pairwise <-  mclapply(group_vals, get_contrast, mc.cores = 16, against = 'Pairwise')

de_results_pairwiseA <- de_results_pairwise %>% bind_rows()


de_table <- bind_rows(de_results) %>% bind_rows(de_results_pairwiseA)
de_table$Organism <-  unique(colData(deseq2_obj)$organism )
write_tsv(de_table, file = args[3])

