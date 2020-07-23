library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

load(args[1]) # umap
load(args[2]) # doublet info

meta <- left_join(umap, doublet_call_table)  %>% 
		mutate(`Doublet Probability` = as.numeric(`Doublet Probability`),
         doublet_score_scran = as.numeric(doublet_score_scran))

subC_remove <- meta %>% 
  group_by(subcluster) %>% 
  summarise(scrublet = mean(`Doublet Probability`), 
            scran = log(mean((doublet_score_scran))), 
            count = n()) %>% 
  arrange(-scran) %>% 
  filter(scrublet > 0.2 | scran > 11, count < 5000) %>% 
  pull(subcluster)

cluster_remove <- meta %>% 
  filter(!subcluster %in% subC_remove) %>% 
  group_by(cluster) %>% 
  summarise(scrublet = mean(`Doublet Probability`), 
            scran = log(mean((doublet_score_scran))), 
            count = n()) %>% 
  arrange(-scran) %>% 
  filter(scrublet > 0.15 | scran > 9, count < 5000) %>% 
  pull(cluster)

cluster_low_n <- meta %>% 
  group_by(cluster) %>% 
  summarise(Count = n()) %>% 
  filter(Count < 20) %>% 
  pull(cluster) 

## for low-n clusters, re-assign to closest cluster in UMAP2D space
## find closest cluster by euclidean dist
#centroid <- meta %>% 
#  group_by(cluster) %>% 
#  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
#dist_obj <- dist(centroid[,2:3]) %>% as.matrix()
#row.names(dist_obj) <- centroid$cluster
#colnames(dist_obj) <- centroid$cluster

#re_label <- list()
#for (i in cluster_low_n){
#  re_label[[i]] <- (dist_obj[,i] %>% sort())[2] %>% names()
#  print((dist_obj[,i] %>% sort())[2] %>% names())
#}

umap <- meta %>% 
  filter(!subcluster %in% subC_remove) %>% 
  filter(!cluster %in% cluster_remove) %>% 
  filter(`Doublet Probability` < 0.8, doublet_score_scran < 1e6) %>% 

save(umap, file = args[3])
