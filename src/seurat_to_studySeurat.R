library(tidyverse)
library(Seurat)
library(SeuratDisk)

load('site/scEiaD_all_seurat_v3.Rdata')

big <- scEiaD

system( 'mkdir -p site/study_level')

#dev_studies <- big@meta.data %>% filter(Tissue == 'Retina') %>% group_by(study_accession, CellType_predict) %>% summarise(Count = n()) %>% mutate(Perc = Count / sum(Count)) %>% filter(Perc > 0.01)  %>% filter(grepl('RPC|Prog|Precur|Neuro', CellType_predict)) %>% pull(study_accession) %>% unique()

for (i in unique(big@meta.data$study_accession)){
#for (i in dev_studies){
	print(i)
	scEiaD <- subset(big, subset = study_accession == i)
	save(scEiaD, file = paste0('site/study_level/', i, '.seurat.Rdata'))
	# h5ad
	## scvelo looks for spliced/unspliced
	scEiaD[['spliced']] <- scEiaD[['RNA']]
	SaveH5Seurat(scEiaD, filename =  paste0('site/study_level/', i, '.h5Seurat'), overwrite = TRUE)
	Convert(paste0('site/study_level/', i, '.h5Seurat'), dest = 'h5ad', overwrite = TRUE)
	system(paste0('rm site/study_level/', i, '.h5Seurat'))
} 

