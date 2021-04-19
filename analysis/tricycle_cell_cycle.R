library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(tricycle)

# tricycle installed with remotes::install_github('hansenlab/tricycle@1ce333ff097661) 
# to avoid the R4.1 requirement
# wget https://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/gene_name_ids.tsv.gz
# wget https://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/scEiaD_all_seurat_v3.Rdata
# the above is 20GB!!!

load('~/data/scEiaD/scEiaD_all_seurat_v3.Rdata')
# scEiaD seurat obj uses the ENSG system for genes names
# so after we turn the seurat object into a sce objecdt
# we swap out the ENGS for the human readble gene names 
# which work better for tricycle
gene_names <- read_tsv('~/data/scEiaD/gene_name_ids.tsv.gz')
sce <- as.SingleCellExperiment(scEiaD_droplet)
row.names(sce) <- row.names(sce) %>% enframe(value = 'ID') %>% left_join(gene_names) %>% pull(Name)

sce <- project_cycle_space(sce, species = 'human', gname.type = 'SYMBOL')

sce <- estimate_cycle_position(sce)

sce <- estimate_cycle_stage(sce, gname.type = 'SYMBOL', species = 'human')

# scVI plot with tricycle phase values
# note the coord flip and y axis flip to match plae.nei.nih.gov
p <- plot_emb_circle_scale(sce, point.size = .1, point.alpha = 0.5, dimred = 'SCVIUMAP') +
  theme_bw(base_size = 14) + 
  coord_flip() + scale_y_reverse()
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
cowplot::plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))

# CCstage (discrete) values plotted on the tricycle PC1/2 embedding projection
scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "CCStage",
                       point_size = 0.2, point_alpha = 0.2) +
  labs(x = "Projected PC1", y = "Projected PC2", title = paste0("Projected cell cycle space (n=", ncol(sce), ")")) +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(size=10))) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 0.9)))

# quick Seurat vs. Tricycle plot
one <- scater::plotReducedDim(sce, dimred = "SCVIUMAP", colour_by = "CCStage",
                       point_size = 0.1, point_alpha = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Tricycle/Schwabe CC Stage") +
  theme_bw(base_size = 14) + 
  coord_flip() + scale_y_reverse() + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 0.9)))

two <- scater::plotReducedDim(sce, dimred = "SCVIUMAP", colour_by = "Phase",
                              point_size = 0.1, point_alpha = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Seurat CC Stage") +
  theme_bw(base_size = 14) + 
  coord_flip() + scale_y_reverse() + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 0.9)))
cowplot::plot_grid(one, two)
