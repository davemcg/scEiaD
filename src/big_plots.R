library(tidyverse)
#library(ggforce)
#library(ggrepel)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

red <- args[1]
load(args[2])

# filter
umap <- umap %>% 
  rename(Stage = integration_group) %>% 
  mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType)) %>% 
  mutate(Size = case_when(organism == 'Homo sapiens' ~ 0.015,
                          TRUE ~ 0.01)) 

# attach colors to cell types
cell_types <- umap %>% 
  pull(CellType) %>% unique() %>% sort()
type_val <- setNames(pals::alphabet(n = cell_types %>% length()), cell_types)
type_col <- scale_colour_manual(values = type_val)
type_fill <- scale_fill_manual(values = type_val)

# cell type known
plot1 <- umap %>% 
  ggplot() + 
  geom_point(aes(x=umap[,paste0(red,'_1')] %>% pull(1), 
                 y = umap[,paste0(red,'_2')] %>% pull(1), 
                 colour = CellType), size = 0.05, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num ), alpha = 0.8, size = 2) +
  type_col + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))

# Age
plot2 <- umap %>% 
  ggplot() + 
  geom_point(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                 y = umap[,paste0(red,'_2')] %>% pull(1),  
                 colour = Stage), size = 0.2, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  scale_color_manual(values = as.vector(pals::alphabet()) %>% tail(2)) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))

# facet by organism
plot3 <- umap %>% 
  ggplot() + 
  geom_point(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                 y = umap[,paste0(red,'_2')] %>% pull(1),   
                 colour = CellType_predict, size = Size), alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_cowplot() + 
  scale_size(guide = 'none') +
  scale_alpha(guide = 'none') +
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  type_col + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~organism) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))


# facet by celltype, color by organism
plot4 <- umap %>% 
  ggplot() + 
  geom_point(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                 y = umap[,paste0(red,'_2')] %>% pull(1), 
                 colour = organism), size = 0.2,  alpha = 0.05) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_cowplot() + 
  scale_size(guide = 'none') +
  scale_alpha(guide = 'none') +
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~CellType) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))

# facet by cluster, color by CellType
plot5 <- umap %>% 
  ggplot() + 
  geom_point(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                 y = umap[,paste0(red,'_2')] %>% pull(1),  
                 colour = CellType), size = 0.2, alpha = 0.05) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_cowplot() + 
  type_col + 
  scale_size(guide = 'none') +
  scale_alpha(guide = 'none') +
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~cluster) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))

png(args[3], width = 1800, height = 5500, res = 150)
plot_grid(plot1, plot4, plot5, ncol = 1, rel_heights = c(0.5,0.5,1))
dev.off()
