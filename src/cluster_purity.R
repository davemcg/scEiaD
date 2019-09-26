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

output_file1 <- args[2]
output_file2 <- args[3]

cluster_purity_plot <- function(obj, algorithm_name){
  obj <- obj %>% filter(!is.na(CellType))
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
  obj <- obj %>% group_by(seurat_clusters, CellType) %>% 
    summarise(Count = n(), CellTypes = list(unique(CellType))) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(seurat_clusters) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) 
  obj$Info <- algorithm_name
  obj <- obj %>% 
    separate(Info, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', 'Dims'), sep = '__')
  list(mean = p1_mean, cluster_count = p1_cluster_count, dist = obj)
}

purity_calcs <- all_obj %>% imap(cluster_purity_plot)


pdf(output_file1, width = 20, height = 4)
#pdf('test.pdf', width = 20, height = 4)
purity_calcs %>% 
  map(3) %>% # third item in list is the purity distributions
  bind_rows() %>% 
  ggplot(aes(x=Method, y = Count, colour = Transform)) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~Covariate + Set + Dims, nrow = 1)
dev.off()

pdf(output_file2, width = 20, height = 4)
#pdf('test2.pdf', width = 20, height = 4)
cbind(purity_calcs %>% 
  map(1) %>% # third item in list is the purity distributions
  enframe(value = 'Mean') %>% 
    unnest(Mean),
  purity_calcs %>% 
    map(2) %>% # third item in list is the purity distributions
    enframe(value = 'Count') %>% unnest(Count) %>% select(Count)) %>% 
  separate(name, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method'), sep = '__') %>% 
  as_tibble() %>% 
  ggplot(aes(x=Method, y = Mean, colour = Transform, size = Count)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~Covariate + Set, nrow = 1)
dev.off()

