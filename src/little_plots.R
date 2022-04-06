library(tidyverse)
library(cowplot)
library(scattermore)
library(glue)
options(rgl.useNULL=TRUE)
args <- commandArgs(trailingOnly = TRUE)
load(args[1])
ptsize = 10
ALPHA=.2

umap <- seurat@meta.data %>% select(-Barcode) %>%
			rename(sceiadUMAP_1 = UMAP_1, sceiadUMAP_2 = UMAP_2) %>% 
			as_tibble(rownames = "Barcode") %>% 
			left_join(seurat@reductions$umap@cell.embeddings %>% as_tibble(rownames = "Barcode"))


plotter <- function(umap,color_by){
	# make labels
	umap$group <- umap[,color_by] %>% pull(1)
	labels <- umap %>% group_by(group) %>%
  			summarise(UMAP_1 = mean(UMAP_1),
           			 UMAP_2 = mean(UMAP_2))
	umap <- umap %>% mutate(group = case_when(is.na(group) ~ 'Missing', TRUE ~ group))
	groups <- umap %>%
	  pull(group) %>% unique() %>% sort()
	type_val <- setNames(c(pals::alphabet(), pals::alphabet2(), pals::glasbey())[1:length(groups)], groups)
	type_col <- scale_colour_manual(values = type_val)
	type_fill <- scale_fill_manual(values = type_val)
	ncells <- nrow(umap)
	study_accession <- umap$study_accession %>% unique()
	sample_accession <- umap$sample_accession %>% unique()
	ggplot_obj <- 	umap %>%
	  ggplot() +
	  geom_scattermore(aes(x = UMAP_1,
                       y = UMAP_2,
                       colour = group), pointsize = (ptsize/5), alpha = ALPHA) +
	  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) +
 	 theme_cowplot() +
	  #geom_label_repel(data = cluster_labels, aes(x=x, y=y, label = seurat_cluster_CellType_num ), alpha = 0.8, size = 2) +
 	 type_col +
	ggrepel::geom_label_repel(data = labels, aes(x=UMAP_1, y=UMAP_2, label = group), alpha = 0.8, size = 2.2) +
 	 theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
 	 xlab('UMAP 1') + ylab(paste('UMAP 2')) +
	  	ggtitle(glue('Study: {study_accession}, Sample: {sample_accession}, Number of Cells: {ncells}'))
	return(ggplot_obj)
}

ct <- plotter(umap, 'CellType')
ctp <- plotter(umap, 'CellType_predict')


png(args[2], width = 2000, height = 3000, res = 200)
plot_grid(ct, ctp,  ncol = 1, rel_heights = c(1,1))
dev.off()
