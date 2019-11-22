library(tidyverse)
library(cowplot)
##########
# load all clustering params
##########
files <- list.files('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/cluster/', 
                    pattern = '*.cluster.Rdata', 
                    full.names = TRUE)
cluster_all <- list()
for (i in files){
  load(i)
  colnames(meta) <- c('Barcode', 'cluster')
  cluster_all[[i]] <- meta
  dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
  knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
  cluster_all[[i]]$dims = dims
  cluster_all[[i]]$knn = knn
}
cluster_all <- cluster_all %>% bind_rows()


########
# load umap file with one cluster param
########
files <- list.files('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/umap/', 
                    pattern = '.*scVI.*mindist0.1.*nneighbors50.*', 
                    full.names = TRUE)
umap_all <- list()
count = 1
for (i in files){
  print(count)
  count = count + 1
  load(i)
  nn = str_extract(i, 'nneighbors\\d+') %>% gsub('nneighbors', '', .) %>% as.numeric()
  dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
  dist = str_extract(i, 'mindist0.\\d+') %>% gsub('mindist', '', .) %>% as.numeric()
  umap_all[[i]] <- umap
  umap_all[[i]]$mindist = dist
  umap_all[[i]]$nneighbors = nn
  umap_all[[i]]$dims = dims
  d <- dims
  for (k in c(cluster_all$knn %>% unique())){
    umap_all[[paste0(i,'_', k)]] <- umap_all[[i]] %>% 
      select(-cluster) %>% 
      left_join(., cluster_all %>% filter(knn == k), by = 'Barcode') %>% 
      mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>%
      group_by(cluster, CellType) %>%
      filter(!is.na(CellType), !is.na(cluster)) %>%
      summarise(Count = n(), x = mean(UMAP_1), y = mean(UMAP_2), 
                Organism = list(unique(organism)), 
                study_accession = list(unique(study_accession)),
                study_n = length(unique(study_accession)),
                org_n = length(unique(Organism)),
                mindist = unique(mindist),
                nneighbors = unique(nneighbors),
                dims = unique(dims),
                Method = unique(Method)) 
    umap_all[[paste0(i,'_', k)]]$knn <- k
  }
  umap_all[[i]] <- NULL
}
umap_one <- umap_all %>% bind_rows() 




# Density Plot of Distribution of Frequency of Top Cell Type in Each Cluster
# Ideally every cluster would be nearly 1 (100%)
umap_one %>% filter(Count > 100) %>% 
  mutate(freq = Count / sum(Count)) %>% 
  ungroup() %>% 
  group_by(dims, cluster, knn) %>% 
  top_n(n = 1, wt = freq) %>% 
  ggplot(aes(freq, group = knn, colour = as.factor(knn))) + geom_density() + theme_cowplot() + facet_wrap(~dims)

# Density Plot of Distribution of Number of Studies in Each Cluster
# Again, higher is better
umap_one %>% filter(Count > 100) %>% 
  mutate(freq = Count / sum(Count)) %>% 
  arrange(cluster, -freq) %>% 
  ungroup() %>% 
  group_by(dims, cluster, knn) %>% 
  top_n(n = 1, wt = freq) %>% 
  ggplot(aes(study_n, group = knn, colour = as.factor(knn))) + geom_density() + theme_cowplot()+ facet_wrap(~dims)


# Density Plot of Distribution of Number of Organisms in Each Cluster
# Again, higher is better
umap_one %>% filter(Count > 100) %>% 
  mutate(freq = Count / sum(Count)) %>% 
  ungroup() %>% 
  group_by(dims, cluster) %>%
  top_n(n = 1, wt = freq) %>% 
  ggplot(aes(org_n, group = dims, colour = as.factor(dims))) + geom_density() + theme_cowplot()
