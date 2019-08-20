args <- commandArgs(trailingOnly = TRUE)

# seurat object
load(args[1])

# metadata
load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')

library(tidyverse)
library(Seurat)
library(schex)
library(cowplot)
library(future)
plan(strategy = "multicore", workers = 12)
options(future.globals.maxSize = 2400000 * 1024^2)
#load('obj.Rdata')


  new_meta <- row.names(obj@meta.data) %>% enframe() %>% mutate(sample_accession = str_extract(value, 'SRS\\d+')) %>% left_join(sra_metadata_extended %>% select(sample_accession, study_accession, Platform, Age, TissueNote) %>% unique()) 
  
  obj@meta.data$Age <- new_meta$Age
  obj@meta.data$Platform <- new_meta$Platform
  obj@meta.data$study_accession <- new_meta$study_accession


###########################
# plot by study
############################
DimPlot(obj, group.by = 'study_accession')
########################
# plot by age
###########################
embeds <- Embeddings(obj[['umap']])
samps <- str_extract(embeds %>% row.names(), '(SRS|iPSC_RPE_scRNA_)\\d+')
new <- samps %>% enframe(value = 'sample_accession') %>% left_join(., sra_metadata_extended %>% select(sample_accession, study_accession, Platform, Age, Source, TissueNote) %>% unique())

embeds <- embeds %>% as_tibble()
embeds <- cbind(embeds, new)
embeds <- embeds %>% mutate(Age = case_when(Age == 1000 ~ 30, TRUE ~ Age))
embeds %>% ggplot(aes(x=UMAP_1, y=UMAP_2, colour = Age)) + 
  geom_point(alpha=0.2, size = 0.1) + 
  scale_color_viridis_c() + theme_minimal()

embeds %>% ggplot(aes(UMAP_1, y=UMAP_2, colour = study_accession)) + 
  geom_point(alpha=0.4, size = 0.1) + 
  scale_color_viridis_d() + 
  theme_minimal() + 
  facet_wrap(~Age)

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
obj <- make_hexbin(obj, nbins = 200)
schex <- list()
for (i in (toupper(c('Rho','Opn1sw', 'Rbpms', 'Sfrp2', 'Olig2', 'Tfap2a', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4')))){
  schex[[i]] <- plot_hexbin_gene(obj, i, action = 'mean', type = 'logcounts') + scale_fill_gradient(low = 'gray', high ='red')
}

cowplot::plot_grid(plotlist = schex)

#FeaturePlot(study_data_integrated, toupper(c('Rho','Opn1sw', 'Rbpms', 'Sfrp2', 'Olig2', 'Tfap2a', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4')), pt.size = 0.1, order= FALSE)