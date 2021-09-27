library(tidyverse)
#library(ggforce)
#library(ggrepel)
library(cowplot)
library(scattermore)
library(glue)
options(rgl.useNULL=TRUE)
args <- commandArgs(trailingOnly = TRUE)
red <- args[1]
load(args[2])
ptsize = 4
ALPHA=.07
if (grepl('onlyWELL', args[2])){
	umapO <- umap
	load(args[4])
	umapO <- umapO %>% left_join(umap %>% select(Barcode, CellType_predict), by = 'Barcode')
	umapO$CellType <- umapO$CellType_predict
	umap <- umapO
 	celltype_col <- 'CellType'
 	print('woo')
    ptsize = 20
 } else { 
 	celltype_col <- 'CellType' 
 }

if (!"cluster" %in% colnames(umap)){
  umap$cluster <- umap$clusters
} 

umapO <- umap
# filter
umap <- umap %>% 
  #rename(Stage = integration_group) %>% 
  mutate(CellType = gsub('Cone Bipolar Cells', 'Bipolar Cells', !!as.symbol(celltype_col))) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets'),
         !grepl('Mesenchyme/Lens', CellType))  %>% 
  mutate(Size = case_when(organism == 'Homo sapiens' ~ 0.015,
                          TRUE ~ 0.01)) 

# attach colors to cell types
cell_types <- umap %>% 
  pull(CellType) %>% unique() %>% sort()
type_val <- setNames(c(pals::alphabet(), pals::alphabet2())[1:length(cell_types)], cell_types)
type_col <- scale_colour_manual(values = type_val)
type_fill <- scale_fill_manual(values = type_val)

# make labels
labels <- umap %>% group_by(CellType) %>% 
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

# cell type known
ncells <- nrow(umap)
plot1 <- umap %>% 
  ggplot() + 
  geom_scattermore(aes(x=umap[,paste0(red,'_1')] %>% pull(1), 
                       y = umap[,paste0(red,'_2')] %>% pull(1), 
                       colour = CellType), pointsize = (ptsize/5), alpha = ALPHA) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num ), alpha = 0.8, size = 2) +
  type_col + 
  ggrepel::geom_label_repel(data = labels, aes(x=UMAP_1, y=UMAP_2, label = CellType ), alpha = 0.8, size = 2.2) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  xlab(paste(red, '1')) + ylab(paste(red, '2')) + 
  ggtitle(glue('Number of Cells: {ncells}'))

# Age
#plot2 <- umap %>% 
#  ggplot() + 
#  geom_scattermore(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
#                       y = umap[,paste0(red,'_2')] %>% pull(1),  
#                       colour = Stage), pointsize = ptsize, alpha = 0.1) + 
#  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
#  theme_cowplot() + 
#  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
#  scale_color_manual(values = as.vector(pals::alphabet()) %>% tail(2)) + 
#  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
#  xlab(paste(red, '1')) + ylab(paste(red, '2'))

# facet by organism
plot3 <- umap %>% 
  ggplot() + 
  geom_scattermore(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                       y = umap[,paste0(red,'_2')] %>% pull(1),   
                       colour = CellType, pointsize = ptsize), alpha = ALPHA) + 
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
  geom_scattermore(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                       y = umap[,paste0(red,'_2')] %>% pull(1), 
                       colour = organism), pointsize = ptsize,  alpha = ALPHA) + 
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
  geom_scattermore(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                       y = umap[,paste0(red,'_2')] %>% pull(1),  
                       colour = CellType), pointsize = ptsize, alpha = ALPHA/2) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_cowplot() + 
  type_col + 
  scale_size(guide = 'none') +
  scale_alpha(guide = 'none') +
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~cluster) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))

# facet by celltype, color by study
plot6 <- umap %>% 
  ggplot() + 
  geom_scattermore(aes(x = umap[,paste0(red,'_1')] %>% pull(1), 
                       y = umap[,paste0(red,'_2')] %>% pull(1), 
                       colour = study_accession), 
                   pointsize = ptsize, alpha = ALPHA) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_cowplot() + 
  scale_size(guide = 'none') +
  scale_alpha(guide = 'none') +
  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num )) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  scale_color_manual(values = pals::alphabet() %>% unname()) +
  facet_wrap(~CellType + organism) +
  xlab(paste(red, '1')) + ylab(paste(red, '2'))

# plot 7
# show both labelled and unlabelled cells, color by study_accession
# purpose: to id studies with "odd" patterning for removal
# e.g. they don't overlap other data
plot7 <- umapO  %>% 
				filter(is.na(TabulaMurisCellType)) %>% 
				mutate(CellType = case_when(is.na(CellType) ~ 'Unknown', TRUE ~ CellType)) %>% 
				mutate(Known = case_when(CellType == 'Unknown' ~ 'Unknown', TRUE ~ 'Known')) %>% 
				ggplot(aes(x=UMAP_1,y=UMAP_2, color = study_accession)) + 
				scattermore::geom_scattermore(pointsize = 0, alpha = 0.8) + 
				scale_color_manual(values = c(pals::alphabet(), pals::alphabet2()) %>% unname()) + cowplot::theme_cowplot() + facet_wrap(~Known) 
png(args[3], width = 5000, height = 18000, res = 300)
plot_grid(plot1, plot4, plot5, plot6, plot7,  ncol = 1, rel_heights = c(0.5,0.5,1, 1,0.5))
dev.off()
