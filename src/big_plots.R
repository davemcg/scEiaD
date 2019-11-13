library(tidyverse)
library(ggforce)
library(ggrepel)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])

# attach colors to cell types
cell_types <- umap %>% mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType)) %>% 
  pull(CellType) %>% unique() %>% sort()
type_val <- setNames(pals::alphabet(n = cell_types %>% length()), cell_types)
type_col <- scale_colour_manual(values = type_val)
type_fill <- scale_fill_manual(values = type_val)


# cell type known
plot1 <- umap %>% 
  mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType)) %>%
  #left_join(., cluster_labels) %>% 
  ggplot() + 
  geom_point(aes(x=UMAP_1, y = UMAP_2, colour = CellType), size = 0.05, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num ), alpha = 0.8, size = 2) +
  type_col + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  xlab('UMAP 1') + ylab('UMAP 2')

# Age
plot2 <- umap %>% 
  rename(Stage = integration_group) %>% 
  mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType)) %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1, y = UMAP_2, colour = Stage), size = 0.01, alpha = 0.01) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  scale_color_manual(values = as.vector(pals::alphabet()) %>% tail(2)) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  xlab('UMAP 1') + ylab('UMAP 2')

# facet by organism
plot3 <- umap %>% 
  mutate(CellType_predict = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>% 
  filter(!is.na(study_accession), !is.na(CellType_predict),
         !CellType_predict %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType_predict)) %>%
  mutate(Size = case_when(organism == 'Homo sapiens' ~ 0.015,
                          TRUE ~ 0.01)) %>% 
  ggplot() + 
  geom_point(aes(x=UMAP_1, y = UMAP_2, colour = CellType_predict, size = Size), alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_cowplot() + 
  scale_size(guide = 'none') +
  scale_alpha(guide = 'none') +
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  type_col + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~organism) +
  xlab('UMAP 1') + ylab('UMAP 2')

pdf(args[2], width = 10, height = 10)
plot_grid(plot1, ncol = 1)
dev.off()
