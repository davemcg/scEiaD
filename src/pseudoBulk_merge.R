args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

output <- args[length(args)]
input <- args[-length(args)]

DE__CELLTYPEPREDICT__res_againstAll <- list()
DE__CELLTYPEPREDICT__res_pairwise <- list()
DE__CELLTYPEPREDICT__res_organism_celltype <- list()
DE__CELLTYPE__res_againstAll <- list()
DE__CELLTYPE__res_pairwise <- list()
DE__CELLTYPE__res_organism_celltype <- list()
DE__CLUSTER__res_againstAll <- list()
DE__CLUSTER__res_pairwise <- list()
DE__CLUSTER__res_organism_celltype <- list()

for (i in input){
	load(i)
	print(i)
	comp <- str_extract(i, '_[ABC][123]_') %>% gsub('_','',.)
	if (comp == 'B1') {
		DE__CELLTYPEPREDICT__res_againstAll[[i]] <- CELLTYPEPREDICT__res_againstAll
	} else if (comp == 'B2') {
		DE__CELLTYPEPREDICT__res_pairwise[[i]] <- CELLTYPEPREDICT__res_pairwise
	} else if (comp == 'B3') {
		DE__CELLTYPEPREDICT__res_organism_celltype[[i]] <- CELLTYPEPREDICT__res_organism_celltype
	} else if (comp == 'A1') {
		DE__CELLTYPE__res_againstAll[[i]] <- CELLTYPE__res_againstAll
	} else if (comp == 'A2') {
		DE__CELLTYPE__res_pairwise[[i]] <- CELLTYPE__res_pairwise
	} else if (comp == 'A3') {
		DE__CELLTYPE__res_organism_celltype[[i]] <- CELLTYPE__res_organism_celltype
	} else if (comp == 'C1') {
		DE__CLUSTER__res_againstAll[[i]] <- CLUSTER__res_againstAll
	} else if (comp == 'C2') {
		DE__CLUSTER__res_pairwise[[i]] <- CLUSTER__res_pairwise
	} else if (comp == 'C3') {
		DE__CLUSTER__res_organism_celltype[[i]] <- CLUSTER__res_organism_celltype
	}
} 

binder <- function(list_obj){
	temp <- list()
	for (i in seq_along(1:length(list_obj))){
		if (sum(grepl('data.frame', class(list_obj[[i]]))) == 1 ){
			temp[[i]] <- list_obj[[i]]
		}
	}
	temp %>% bind_rows()
}

PB_resultsABC <- list()
PB_resultsC2 <- list()
if (grepl('ABC', output)) {
	DE__CELLTYPEPREDICT__res_againstAll <- binder(DE__CELLTYPEPREDICT__res_againstAll) %>% mutate(PB_Test = 'CellType (Predict) against Remaining')
	DE__CELLTYPEPREDICT__res_pairwise <- binder(DE__CELLTYPEPREDICT__res_pairwise) %>% mutate(PB_Test = 'Pairwise CellType (Predict) against CellType (Predict)')
	DE__CELLTYPEPREDICT__res_organism_celltype <- binder(DE__CELLTYPEPREDICT__res_organism_celltype) %>% mutate(PB_Test = 'Organism against Organism within CellType (Predict)')
	DE__CELLTYPE__res_againstAll <- binder(DE__CELLTYPE__res_againstAll) %>% mutate(PB_Test = 'CellType against Remaining') 
	DE__CELLTYPE__res_pairwise <- binder(DE__CELLTYPE__res_pairwise) %>% mutate(PB_Test = 'Pairwise CellType against CellType')
	DE__CELLTYPE__res_organism_celltype <- binder(DE__CELLTYPE__res_organism_celltype) %>% mutate(PB_Test = 'Organism against Organism within CellType')
	DE__CLUSTER__res_againstAll <- binder(DE__CLUSTER__res_againstAll) %>% mutate(PB_Test = 'Cluster against Remaining')
	DE__CLUSTER__res_organism_celltype <- binder(DE__CLUSTER__res_organism_celltype) %>% mutate(PB_Test = 'Organism against Organism within Cluster')

	PB_resultsABC <- bind_rows(DE__CELLTYPEPREDICT__res_againstAll,
									DE__CELLTYPEPREDICT__res_pairwise,
									DE__CELLTYPEPREDICT__res_organism_celltype,
									DE__CELLTYPE__res_againstAll,
									DE__CELLTYPE__res_pairwise,
									DE__CELLTYPE__res_organism_celltype,
									DE__CLUSTER__res_againstAll,
									DE__CLUSTER__res_organism_celltype)
	save(PB_resultsABC, file = output)
} else {
	DE__CLUSTER__res_pairwise <- binder(DE__CLUSTER__res_pairwise) %>% mutate(PB_Test = 'Pairwise Cluster against Cluster')
	PB_resultsC2 <- DE__CLUSTER__res_pairwise
	save(PB_resultsC2, file = output)
}

