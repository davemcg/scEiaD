args <- commandArgs(trailingOnly = TRUE)

# umap objects
load(args[1])

library(tidyverse)
library(Polychrome)
library(schex)
library(cowplot)




###########################
# plot by study
############################

naive  %>% ggplot(aes(x=UMAP_1, y=UMAP_2, color = study_accession)) + 
  geom_point(size = 0.01, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  ggsci::scale_color_npg() + ggtitle('No Integration')

p1 <- seurat %>% ggplot(aes(x=UMAP_1, y=UMAP_2, color = study_accession)) + 
  geom_point(size = 0.01, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  ggsci::scale_color_npg() + ggtitle('Seurat CCA')

p2 <- harmony %>% ggplot(aes(x=UMAP_1, y=UMAP_2, color = study_accession)) + 
  geom_point(size = 0.01, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  ggsci::scale_color_npg() + ggtitle('Harmony')


p3 <- mnn %>% ggplot(aes(x=UMAP_1, y=UMAP_2, color = study_accession)) + 
  geom_point(size = 0.01, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  ggsci::scale_color_npg() + ggtitle('MNN')

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 20))
)
plot_grid(p1 + theme(legend.position="none"), 
          p2 + theme(legend.position="none"), 
          p3 + theme(legend.position="none"), 
          legend,
          rel_widths = c(2,2,2,0.6),
          nrow = 1)


########################
# plot facetted by age, colored by study
###########################
p1 <- seurat %>% ggplot(aes(UMAP_1, y=UMAP_2, colour = study_accession)) + 
  geom_point(alpha=0.4, size = 0.1) + 
  scale_color_viridis_d() + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_minimal() + 
  facet_wrap(~Age)
p1

p2 <- harmony %>% ggplot(aes(UMAP_1, y=UMAP_2, colour = study_accession)) + 
  geom_point(alpha=0.4, size = 0.1) + 
  scale_color_viridis_d() + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_minimal() + 
  facet_wrap(~Age)
p2

p3 <- mnn %>% ggplot(aes(UMAP_1, y=UMAP_2, colour = study_accession)) + 
  geom_point(alpha=0.4, size = 0.1) + 
  scale_color_viridis_d() + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_minimal() + 
  facet_wrap(~Age)
p3 

#########################
# plot facetted by age, colored by inferred ID
##########################
p1 <- seurat %>% filter(!new_CellType_transfer %in% c('Doublets', 'Red Blood Cells')) %>% 
  ggplot(aes(UMAP_1, y=UMAP_2, colour = new_CellType_transfer)) + 
  geom_point(alpha=0.4, size = 0.1) +  theme_black() + 
  scale_color_manual(values = unname(alphabet.colors())) +
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~Age)
p1

p2 <- harmony %>% filter(!new_CellType_transfer %in% c('Doublets', 'Red Blood Cells')) %>% 
  ggplot(aes(UMAP_1, y=UMAP_2, colour = new_CellType_transfer)) + 
  geom_point(alpha=0.4, size = 0.1) +  theme_black() + 
  scale_color_manual(values = unname(alphabet.colors())) +
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~Age)
p2

p3 <- mnn %>% filter(!new_CellType_transfer %in% c('Doublets', 'Red Blood Cells')) %>% 
  ggplot(aes(UMAP_1, y=UMAP_2, colour = new_CellType_transfer)) + 
  geom_point(alpha=0.4, size = 0.1) +  theme_black() + 
  scale_color_manual(values = unname(alphabet.colors())) +
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~Age)
p3

###################################
# plot facetted by seurat_cluster, colored by inferred ID
####################################
p1 <- seurat %>% filter(!new_CellType_transfer %in% c('Doublets', 'Red Blood Cells')) %>% 
  mutate(integrated_snn_res.1 = as.numeric(as.character(integrated_snn_res.1))) %>% 
  ggplot(aes(UMAP_1, y=UMAP_2, colour = new_CellType_transfer)) + 
  geom_point(alpha=0.4, size = 0.1) +  theme_black() + 
  scale_color_manual(values = unname(alphabet.colors())) +
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~integrated_snn_res.1)
p1

p2 <- harmony %>% filter(!new_CellType_transfer %in% c('Doublets', 'Red Blood Cells')) %>% 
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters))) %>% 
  ggplot(aes(UMAP_1, y=UMAP_2, colour = new_CellType_transfer)) + 
  geom_point(alpha=0.4, size = 0.1) +  theme_black() + 
  scale_color_manual(values = unname(alphabet.colors())) +
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~seurat_clusters)
p2

p3 <- mnn %>% filter(!new_CellType_transfer %in% c('Doublets', 'Red Blood Cells')) %>% 
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters))) %>% 
  ggplot(aes(UMAP_1, y=UMAP_2, colour = new_CellType_transfer)) + 
  geom_point(alpha=0.4, size = 0.1) +  theme_black() + 
  scale_color_manual(values = unname(alphabet.colors())) +
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~seurat_clusters)
p3


#########################################
# rough assessment of cluster purity
##########################################
seurat %>%
  mutate(integrated_snn_res.1 = as.numeric(as.character(integrated_snn_res.1))) %>% 
  group_by(integrated_snn_res.1, new_CellType_transfer) %>% summarise(Count = n()) %>% mutate(Perc = Count/sum(Count)) %>% filter(Perc > 0.1)

harmony %>% 
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters))) %>% 
  group_by(seurat_clusters, new_CellType_transfer) %>% summarise(Count = n()) %>% mutate(Perc = Count/sum(Count)) %>% filter(Perc > 0.1)

mnn %>% 
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters))) %>% 
  group_by(seurat_clusters, new_CellType_transfer) %>% summarise(Count = n()) %>% mutate(Perc = Count/sum(Count)) %>% filter(Perc > 0.1)

p1_mean <- seurat %>% group_by(integrated_snn_res.1, new_CellType_transfer) %>% 
  summarise(Count = n()) %>% 
  mutate(Perc = Count/sum(Count)) %>% 
  filter(Perc > 0.1) %>% 
  ungroup() %>% group_by(integrated_snn_res.1) %>% 
  summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
  pull(Count) %>% mean() %>% round(., 3)
p1 <- seurat %>% group_by(integrated_snn_res.1, new_CellType_transfer) %>% 
  summarise(Count = n()) %>% 
  mutate(Perc = Count/sum(Count)) %>% 
  filter(Perc > 0.1) %>% 
  ungroup() %>% group_by(integrated_snn_res.1) %>% 
  summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
  ggplot(aes(x=paste('Seurat, mean =', p1_mean), y = Count)) +
  geom_violin() +
  ggbeeswarm::geom_quasirandom() + xlab('') 


p2_mean <- harmony %>% group_by(seurat_clusters, new_CellType_transfer) %>% 
  summarise(Count = n()) %>% 
  mutate(Perc = Count/sum(Count)) %>% 
  filter(Perc > 0.1) %>% 
  ungroup() %>% group_by(seurat_clusters) %>% 
  summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
  pull(Count) %>% mean() %>% round(., 3)
p2 <- harmony %>% group_by(seurat_clusters, new_CellType_transfer) %>% 
  summarise(Count = n()) %>% 
  mutate(Perc = Count/sum(Count)) %>% 
  filter(Perc > 0.1) %>% 
  ungroup() %>% group_by(seurat_clusters) %>% 
  summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
  ggplot(aes(x=paste('Harmony, mean =', p2_mean), y = Count)) +
  geom_violin() +
  ggbeeswarm::geom_quasirandom() + xlab('') 

p3_mean <- mnn %>% group_by(seurat_clusters, new_CellType_transfer) %>% 
  summarise(Count = n()) %>% 
  mutate(Perc = Count/sum(Count)) %>% 
  filter(Perc > 0.1) %>% 
  ungroup() %>% group_by(seurat_clusters) %>% 
  summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
  pull(Count) %>% mean() %>% round(., 3)
p3 <- mnn %>% group_by(seurat_clusters, new_CellType_transfer) %>% 
  summarise(Count = n()) %>% 
  mutate(Perc = Count/sum(Count)) %>% 
  filter(Perc > 0.1) %>% 
  ungroup() %>% group_by(seurat_clusters) %>% 
  summarise(Count = n()) %>% 
  ggplot(aes(x=paste('MNN, mean =', p3_mean), y = Count)) +
  geom_violin() +
  ggbeeswarm::geom_quasirandom() + xlab('')

plot_grid(p1, p2, p3, nrow = 1)  

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
markers <- list()
for (i in (toupper(c('Rho','Opn1sw', 'Rbpms', 'Sfrp2', 'Olig2', 'Tfap2a', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4', 'Dct')))){
  markers[[i]] <- plot_hexbin_gene(obj, i, action = 'mean', type = 'logcounts') + scale_fill_gradient(low = 'gray', high ='blue')
}

cowplot::plot_grid(plotlist = markers)

#FeaturePlot(study_data_integrated, toupper(c('Rho','Opn1sw', 'Rbpms', 'Sfrp2', 'Olig2', 'Tfap2a', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4')), pt.size = 0.1, order= FALSE)

############ post integration 
# cluster purity
meta %>% as_tibble(rownames = 'Barcode') %>% left_join(.,anno, by = 'Barcode')  %>% 
  ggplot(aes(x=`UMAP_1.x`, y = `UMAP_2.x`, colour = `ID`)) + 
  geom_point(size=0.01) + 
  theme_minimal()  + 
  guides(colour = guide_legend(override.aes = list(size = 2)))

meta %>% as_tibble(rownames = 'Barcode') %>% left_join(.,anno, by = 'Barcode')  %>% 
  ggplot(aes(x=`UMAP_1.x`, y = `UMAP_2.x`, colour = as.factor(`Age.x`))) + 
  geom_point(size=0.01) + scale_color_viridis_d() +
  theme_minimal()  + facet_wrap(~ID) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

