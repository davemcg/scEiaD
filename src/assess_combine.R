args <- commandArgs(trailingOnly = TRUE)

# seurat object
load(args[1])

# metadata
load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')

library(tidyverse)
library(Seurat)

###########################
# plot by study
############################

########################
# plot by age
###########################

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