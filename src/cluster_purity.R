# cluster purity plot

# the idea is that, theoretically..., every seurat cluster is one cell type
# so let's quantify that with the known labels
args <- commandArgs(trailingOnly = TRUE)

library(cowplot)
library(tidyverse)

umap_folder <- args[1]
umap_files <- list.files(umap_folder, "*umap.Rdata", full.names = TRUE)
all_obj <- list()
for (i in umap_files){
  load(i)
  set <- gsub('umap/|.umap.Rdata','', i)
  all_obj[[set]] <- umap
}

output_file <- args[2]


cluster_purity_plot <- function(obj, algorithm_name){
  p1_mean <- obj %>% group_by(seurat_clusters, CellType) %>% 
    summarise(Count = n()) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(seurat_clusters) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
    pull(Count) %>% mean() %>% round(., 2)
  p1_cluster_count <- obj %>% group_by(seurat_clusters, CellType) %>% 
    summarise(Count = n()) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(seurat_clusters) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
    pull(seurat_clusters) %>% unique() %>% length()
  obj %>% group_by(seurat_clusters, CellType) %>% 
    summarise(Count = n()) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(seurat_clusters) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
    ggplot(aes(x=paste0(algorithm_name,' mean = ', p1_mean, ', Sum = ', p1_cluster_count), y = Count)) +
    geom_violin() +
    ggbeeswarm::geom_quasirandom() + xlab('') +
    coord_cartesian(ylim = c(1,5))
}

purity_plots <- list()
x <- all_obj %>% map(cluster_purity_plot, algorithm_name = 'bloo')

p1 <- cluster_purity_plot(seurat,'Seurat')
p2 <- cluster_purity_plot(harmony,'Harmony')
p3 <- cluster_purity_plot(mnn,'MNN')
p4 <- cluster_purity_plot(naive,'Naive')

plot_grid(p4, p1, p2, p3, nrow = 1)  