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

# load counts
load('seurat_obj/raw/n_features-500__transform-counts__partition-PR__covariate-batch__preFilter.seuratV3.Rdata')

# opn1lw == M/Lcones
mlCones <- seurat__standard@assays$RNA@counts['ENSG00000102076',] %>% enframe(name = 'Barcode') %>% filter(value > 1)
# opn1sw == Scones
sCones <- seurat__standard@assays$RNA@counts['ENSG00000128617',] %>% enframe(name = 'Barcode') %>% filter(value > 1)
# rho == rods
Rods <- seurat__standard@assays$RNA@counts['ENSG00000163914',] %>% enframe(name = 'Barcode') %>% filter(value > 1)
# arr3 == cones
Cones <- seurat__standard@assays$RNA@counts['ENSG00000120500',] %>% enframe(name = 'Barcode') %>% filter(value > 1)

plotter <- function(file){
	load(file)
	umap <- umap %>%
		mutate(PR = case_when(
							Barcode %in% mlCones$Barcode ~ 'OPN1LW',
							Barcode %in% sCones$Barcode ~ 'OPN1SW',
							Barcode %in% Rods$Barcode ~ 'RHO',
							Barcode %in% Cones$Barcode ~ 'ARR3'))
	plot1 <- umap %>% 
 	 ggplot() + 
	  geom_scattermore(aes(x=umap[,paste0(red,'_1')] %>% pull(1), 
                       y = umap[,paste0(red,'_2')] %>% pull(1), 
                       colour = PR), pointsize = 4, alpha = ALPHA) + 
 	 guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
	  theme_cowplot() +
	  ggsci::scale_color_aaas() + 
	  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
	  xlab(paste(red, '1')) + ylab(paste(red, '2')) + 
	  ggtitle(file) + facet_wrap(~cluster)
	plot1
}

cluster_PR_type <- function(file){
	load(file)
	umap <- umap %>%
		mutate(PR = case_when(
							Barcode %in% mlCones$Barcode ~ 'OPN1LW',
							Barcode %in% sCones$Barcode ~ 'OPN1SW',
							Barcode %in% Rods$Barcode ~ 'RHO'
							))
	opn1sw_good_cluster_count <- umap %>% filter(!is.na(PR)) %>% group_by(cluster, PR) %>% summarise(Count = n()) %>% mutate(Ratio = Count / sum(Count)) %>% filter(PR == 'OPN1SW', Ratio > 0.7, Count > 500) %>% nrow()
	if (opn1sw_good_cluster_count > 0){
		print(file)
		print(umap %>% filter(!is.na(PR)) %>% group_by(cluster, PR) %>% summarise(Count = n()) %>% mutate(Ratio = Count / sum(Count)) %>% filter(Ratio > 0.5) %>% data.frame())
		print(umap %>% filter(!is.na(PR)) %>% group_by(cluster, PR) %>% summarise(Count = n()) %>% mutate(Ratio = Count / sum(Count)) %>% filter(Ratio > 0.5) %>% mutate(tot = Count * Ratio) %>% pull(tot) %>%  sum())
	}
}

for (i in files){
	print(i)
	name <- gsub('__umap.Rdata','', basename(i))
	name <- glue::glue('pipeline_data/plots/PR/{name}.pdf')
	pdf(name)
	print(plotter(i))
	dev.off()
}

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
