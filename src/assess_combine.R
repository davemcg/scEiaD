args <- commandArgs(trailingOnly = TRUE)

# seurat object
load(args[1])

# metadata
load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')

library(tidyverse)
library(Seurat)
load('seurat_merged.Rdata')
new_meta <- row.names(seurat_merged@meta.data) %>% enframe() %>% mutate(sample_accession = str_extract(value, 'SRS\\d+')) %>% left_join(sra_metadata_extended %>% select(sample_accession, study_accession, Platform, Age, TissueNote) %>% unique()) 

seurat_merged@meta.data$Age <- new_meta$Age
seurat_merged@meta.data$Platform <- new_meta$Platform
seurat_merged@meta.data$study_accession <- new_meta$study_accession

###########################
# plot by study
############################

########################
# plot by age
###########################


new_meta <- row.names(seurat_merged@meta.data) %>% enframe() %>% mutate(sample_accession = str_extract(value, 'SRS\\d+')) %>% left_join(sra_metadata_extended %>% select(sample_accession, study_accession, Platform, Age, TissueNote) %>% unique()) 

seurat_merged@meta.data$Age <- new_meta$Age
seurat_merged@meta.data$Platform <- new_meta$Platform
seurat_merged@meta.data$study_accession <- new_meta$study_accession


########################
# plot by markers
# Rho for rods
# Opn1sw for cones
# Sfrp2 for early progenitors
# Olig2 for neurogenic
# Tfap2a for amacrine
# Ccnd1 for late progenitors
# Aqp4 for muller glia
# Vsx1 to bipolar
# Elavl4 for RGC
####################################
FeaturePlot(study_data_integrated, toupper(c('Rho','Opn1sw', 'Rbpms', 'Sfrp2', 'Olig2', 'Tfap2a', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4')), pt.size = 0.1, order= FALSE)