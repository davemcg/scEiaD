library(tidyverse)
library(cowplot)

# load annotations
# load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/cell_info_labelled.Rdata')

x2 = '.*n_features2000'
x5 = '.*n_features5000'
x10 = '.*n_features10000'
umap_mega <- list()
for (w in c('scVI','fast')){
  for (x in c(x2,x5, x10)){
    # load all clustering params -----
    files <- list.files('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/cluster/', 
                        pattern = paste0(x, '.*', method, '.*.cluster.Rdata'), 
                        full.names = TRUE)
    cluster_all <- list()
    for (i in files){
      load(i)
      colnames(meta) <- c('Barcode', 'cluster')
      cluster_all[[i]] <- meta
      nf = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
      dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
      knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
      method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
      cluster_all[[i]]$dims = dims
      cluster_all[[i]]$knn = knn
      cluster_all[[i]]$nf = nf
      cluster_all[[i]]$cluster_sum = max(meta$cluster)
      cluster_all[[i]]$median_cluster_n = max(meta$cluster)
      cluster_all[[i]]$method <- method
    }
    cluster_all <- cluster_all %>% bind_rows()
    
    
    
    # load umap file with one cluster param -----
    files <- list.files('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/umap/', 
                        pattern = paste0(x,'.*scVI.*mindist0.1.*nneighbors50.*'), 
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
      method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
      umap_all[[i]] <- umap
      umap_all[[i]]$mindist = dist
      umap_all[[i]]$nneighbors = nn
      umap_all[[i]]$dims = dims
      umap_all[[i]]$method = method
      d <- dims
      for (k in c(cluster_all$knn %>% unique())){
        umap_all[[paste0(i,'_', k)]] <- umap_all[[i]] %>% 
          select(-cluster) %>% 
          left_join(., cluster_all %>% filter(dims == d, knn == k), by = 'Barcode') %>% 
          filter(!is.na(CellType), !is.na(cluster)) %>%
          mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>%
          group_by(cluster, CellType, organism, study_accession, method) %>% 
          summarise(Count = n(),
                    x = mean(UMAP_1), 
                    y = mean(UMAP_2),
                    mindist = unique(mindist),
                    nneighbors = unique(nneighbors),
                    dims = unique(dims),
                    Method = unique(Method)) %>% 
          summarise(study_count = n(), 
                    Count = sum(Count),
                    x = mean(x), 
                    y = mean(y),
                    mindist = unique(mindist),
                    nneighbors = unique(nneighbors),
                    dims = unique(dims),
                    Method = unique(Method)) %>% 
          summarise(organism_count = n(), 
                    Count = sum(Count), 
                    study_count = sum(study_count),
                    x = mean(x), 
                    y = mean(y),
                    mindist = unique(mindist),
                    nneighbors = unique(nneighbors),
                    dims = unique(dims),
                    Method = unique(Method)) %>% 
          mutate(freq = Count / sum(Count)) %>% 
          ungroup() %>% 
          group_by(dims, cluster) %>% 
          top_n(n = 1, wt = freq) 
        umap_all[[paste0(i,'_', k)]]$knn <- k
        
      }
      umap_all[[i]] <- NULL
    }
    umap_one <- umap_all %>% bind_rows() 
    umap_one$nfeatures <- x
    umap_mega[[x]] <- umap_one
  }
}
umap_one <- umap_mega %>% bind_rows()

# boxplot Plot of Distribution of Frequency of Top Cell Type in Each Cluster
# Ideally every cluster would be nearly 1 (100%)
umap_one %>% filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric()) %>% 
  ggplot(aes(y = freq, x = as.factor(dims), fill = as.factor(knn))) + 
  geom_boxplot() + geom_hline(aes(yintercept=c(0.855))) +
  theme_cowplot() + facet_wrap(~nfeatures, ncol = 1, strip.position="right") + coord_flip()

## pointrange
p1 <- umap_one %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% group_by(knn, nfeatures, dims) %>% 
  summarise(mean = mean(freq), 
            median = median(freq),
            q1 = quantile(freq, 0.25),
            q3 = quantile(freq, 0.75)) %>% 
  ggplot(aes(y = mean, x = as.factor(dims), 
             color = as.factor(knn))) + 
  geom_pointrange(aes(y=mean, ymin = q1, ymax = q3), position = position_dodge(width = 0.5), show.legend = FALSE) + 
  theme_cowplot() + coord_flip() + scale_color_brewer(palette = 'Set2') +
  facet_wrap(~nfeatures, ncol = 1, strip.position="right") +
  ylab('Mean Cell Type Purity Per Cluster')


# Density Plot of Distribution of Number of Studies in Each Cluster
# Again, higher is better
umap_one %>% filter(Count > 1000) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric()) %>% 
  ggplot(aes(y = study_count, x = as.factor(dims), fill = as.factor(knn))) + 
  geom_boxplot() +
  theme_cowplot() + facet_wrap(~nfeatures, ncol = 1, strip.position="right") + coord_flip()
## pointrange
p2 <- umap_one %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% group_by(knn, nfeatures, dims) %>% 
  summarise(mean = mean(study_count), 
            median = median(study_count),
            q1 = quantile(study_count, 0.25),
            q3 = quantile(study_count, 0.75)) %>% 
  ggplot(aes(y = mean, x = as.factor(dims), 
             color = as.factor(knn))) + 
  geom_pointrange(aes(y=mean, ymin = q1, ymax = q3), 
                  position = position_dodge(width = 0.5), show.legend = FALSE) + 
  theme_cowplot() + coord_flip() + scale_color_brewer(palette = 'Set2') +
  facet_wrap(~nfeatures, ncol = 1, strip.position="right") +
  ylab('Mean Number of Studies\nPer Cluster/Cell Type')


# Density Plot of Distribution of Number of Organisms in Each Cluster
# Again, higher is better
umap_one %>% filter(Count > 1000) %>% 
  #mutate(org_n = length(Organism[[1]])) %>% 
  ggplot(aes(organism_count, group = knn, colour = as.factor(knn))) + geom_density() + theme_cowplot() + facet_wrap(~dims+nfeatures)

## pointrange
p3 <- umap_one %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% group_by(knn, nfeatures, dims) %>% 
  summarise(mean = mean(organism_count), 
            median = median(organism_count),
            q1 = quantile(organism_count, 0.25),
            q3 = quantile(organism_count, 0.75)) %>% 
  ggplot(aes(y = mean, x = as.factor(dims), 
             color = as.factor(knn))) + 
  geom_pointrange(aes(y=mean, ymin = q1, ymax = q3), position = position_dodge(width = 0.5)) + 
  theme_cowplot() + coord_flip() + scale_color_brewer(palette = 'Set2') +
  facet_wrap(~nfeatures, ncol = 1, strip.position="right") +
  ylab('Mean Number of Organisms\nPer Cluster/Cell Type')


patchwork::plot_layout(p1 + p2 + p3, nrow = 1)


# break down metrics by most represented (by study n) cell types
p4 <- umap_one %>% 
  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods', 
                         'Cones', 'Amacrine Cells', 'Microglia', 
                         'Pericytes')) %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% 
  group_by(knn, nfeatures, dims, CellType) %>% 
  summarise(mean = mean(freq), 
            median = median(freq)) %>% 
  ggplot(aes(y = mean, x = as.factor(dims), color = as.factor(knn)), group = CellType) + 
  geom_point(show.legend = FALSE) +
  theme_cowplot() + 
  facet_grid(rows = vars(CellType), cols = vars(nfeatures)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p5 <- umap_one %>% filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods', 
                                          'Cones', 'Amacrine Cells', 'Microglia', 
                                          'Pericytes')) %>% 
  group_by(CellType, dims, knn, nfeatures, study_count) %>% 
  summarise(Count = sum(Count)) %>%  
  mutate(freq = Count / sum(Count)) %>% 
  summarise(study_mean = sum(study_count * freq)) %>%  
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric()) %>% 
  ggplot(aes(y = study_mean, x = as.factor(dims), color = as.factor(knn))) + 
  geom_point(show.legend = FALSE) +
  theme_cowplot() + 
  facet_grid(rows = vars(CellType), cols = vars(nfeatures)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# average number of organisms 
# present in a cluster for each cell type
p6 <- umap_one %>% filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods', 
                                          'Cones', 'Amacrine Cells', 'Microglia', 
                                          'Pericytes')) %>% 
  group_by(CellType, dims, knn, nfeatures, organism_count) %>% 
  summarise(Count = sum(Count)) %>%  
  mutate(freq = Count / sum(Count)) %>% 
  summarise(organism_mean = sum(organism_count * freq)) %>%  
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric()) %>% 
  ggplot(aes(y = organism_mean, x = as.factor(dims), color = as.factor(knn))) + 
  geom_point() +
  theme_cowplot() + 
  facet_grid(rows = vars(CellType), cols = vars(nfeatures)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

patchwork::plot_layout(p4 + p5 + p6, nrow = 1)





# sum sore metric (averaged across all labelled cells)
celltype <- umap_one %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% group_by(knn, nfeatures, dims) %>% 
  summarise(celltype_mean = mean(freq), 
            celltype_median = median(freq))
study <- umap_one %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% group_by(knn, nfeatures, dims) %>% 
  summarise(study_mean = mean(study_count), 
            study_median = median(study_count))
org <-  umap_one %>% 
  filter(Count > 100) %>% 
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% 
           as.numeric()) %>% group_by(knn, nfeatures, dims) %>% 
  summarise(org_mean = mean(organism_count), 
            org_median = median(organism_count))

tab_overall <- left_join(celltype, org) %>% 
  left_join(study) %>% 
  select(-contains("median")) %>% 
  ungroup() %>% 
  mutate(Sum_mean = #(org_mean/max(org_mean)) + 
           (celltype_mean/max(celltype_mean)) + 
           (study_mean/max(study_mean)),
         org = org_mean/max(org_mean),
         ct = celltype_mean/max(celltype_mean),
         sm = study_mean/max(study_mean)) %>% 
  arrange(-Sum_mean)




# sum score metric (scaled by highest mean in each category)
# broken down by represented cell type
# each cell type equally represented 
org_mean <- umap_one %>% filter(Count > 100) %>% 
  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods', 
                         'Cones', 'Amacrine Cells', 'Microglia', 
                         'Pericytes')) %>% 
  group_by(CellType, dims, knn, nfeatures, organism_count) %>% 
  summarise(Count = sum(Count)) %>%  
  mutate(freq = Count / sum(Count)) %>% 
  summarise(organism_mean = (sum(organism_count * freq)) / 3) %>%  
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric())

study_mean <- umap_one %>% filter(Count > 100) %>% 
  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods', 
                         'Cones', 'Amacrine Cells', 'Microglia', 
                         'Pericytes')) %>% 
  group_by(CellType, dims, knn, nfeatures, study_count) %>% 
  summarise(Count = sum(Count)) %>%  
  mutate(freq = Count / sum(Count)) %>% 
  summarise(study_mean = (sum(study_count * freq)) / 9) %>%  
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric())

celltype_mean <- umap_one %>% filter(Count > 100) %>% 
  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods', 
                         'Cones', 'Amacrine Cells', 'Microglia', 
                         'Pericytes')) %>% 
  group_by(CellType, dims, knn, nfeatures, CellType) %>% 
  summarise(Count = sum(Count)) %>%  
  mutate(freq = Count / sum(Count)) %>% 
  # mutate(celltype_mean = sum(Count * freq) / sum(Count)) %>%  
  mutate(nfeatures = str_extract(nfeatures, '\\d+') %>% as.numeric())

tab_celltype <- left_join(org_mean, celltype_mean) %>% 
  left_join(study_mean) %>% 
  mutate(Sum = organism_mean + freq + study_mean) %>% 
  ungroup() %>% 
  group_by(dims, knn, nfeatures) %>% 
  summarise(Sum = sum(Sum),
            organism = sum(organism_mean),
            freq = sum(freq),
            study = sum(study_mean)) %>% 
  ungroup() %>% 
  mutate(organism = organism / max(organism),
         freq = freq / max(freq),
         study = study/ max(study)) %>% 
  mutate(Sum = freq + study)

left_join(tab_overall, tab_celltype) %>% mutate(SumSum = Sum_mean + Sum) %>% arrange(-SumSum) %>% filter(dims != 8)
