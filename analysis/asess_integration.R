# celltype prediction accuracy
load('~/data/scEiaD_2022_02/results/merged_stats.Rdata')
load('~/data/scEiaD_2022_02/results/merged_xgboost.Rdata')
library(tidyverse)
library(cowplot)
accuracy <- accuracy %>% 
  filter(score != 0) %>% 
  mutate(
    partition = str_extract(file, 'partition-[a-zA-Z]+') %>% gsub('partition-','',.),
    norm = str_extract(file, 'transform-[a-zA-Z]+') %>% gsub('transform-','',.),
    nf = str_extract(file, 'n_features-\\d+') %>% gsub('n_features-', '', .) %>% as.numeric(),
    dims = str_extract(file, 'dims-\\d+') %>% gsub('dims-', '', .) %>% as.numeric(),
    method = str_extract(file, 'method-[a-zA-Z]+') %>% gsub('method-','',.),
    knn = str_extract(file, 'knn-\\d+\\.\\d+|knn-\\d+') %>% gsub('knn-', '', .) %>% as.numeric(),
    epochs = str_extract(file, 'epochs-\\d+') %>% gsub('epochs-','',.) %>% as.numeric())

purity <- 1; balance <- 1; mixing <- 1
partition = 'universe'
method = 'scvIprojection'

accuracy %>% 
  filter(partition %in% c('universe')) %>% 
  group_by(norm, nf, dims, method, knn, epochs) %>% 
  summarise(bad_score_sum = sum(score < 0.8), xgboost_score = mean(score)) %>% 
  arrange(-xgboost_score)

perf %>% 
  filter(partition %in% c('universe')) %>% 
  select(-clusterN, -clusterMedian, -subset) %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  select(-LISI_CellType) %>% 
  left_join(accuracy %>% 
              filter(partition %in% c('universe')) %>% 
              group_by(norm, nf, dims, method, knn, epochs) %>% 
              summarise(bad_score_sum = sum(score < 0.8), xgboost_score = mean(score))) %>% 
  left_join(perf_CT %>% 
              bind_rows() %>% 
              filter(Score == 'Silhouette') %>% 
              group_by(normalization, nf, dims, method, knn, epochs) %>% 
              summarise(Silhouette_CellType = mean(Value, na.rm = TRUE)) 
              ) %>% 
  left_join(perf_CT %>% 
              bind_rows() %>% 
              filter(Score == 'LISI') %>% 
              group_by(normalization, nf, dims, method, knn, epochs) %>% 
              summarise(LISI_CellType = mean(Value, na.rm = TRUE)) 
  ) %>% 
  mutate(Score =
           scales::rescale(LISI_CellType, c(1,0)) * purity +
           scales::rescale(Silhouette_CellType, c(0,1)) * purity  +
           scales::rescale(Silhouette_Cluster, c(0,1)) * purity  +
           scales::rescale(LISI_Cluster, c(1,0)) * purity  +
           
           scales::rescale(`xgboost_score`, c(0,1)) * balance +
           scales::rescale(`NMI_CellType-Cluster`, c(0,1)) * balance  +
           scales::rescale(`ARI_CellType-Cluster`, c(0,1)) * balance  +
           
           scales::rescale(LISI_Batch, c(0,1)) * mixing  + # Z score
           scales::rescale(Silhouette_Batch, c(1,0)) * mixing,
         Score = Score / 9) %>% 
  mutate(nf = as.factor(nf),
         `scVI latent dims` = as.factor(dims),
         method = gsub('scVIprojection', 'scVI-projection',method),
         method = gsub('scVI$', 'scVI-standard', method),
         method = factor(method, levels = c('scVI-standard','scVI-projection'))) %>% 
  
ggplot(aes(x=Score, y = nf, color = `scVI latent dims`, shape = as.factor(epochs))) + 
  facet_wrap(~Run) +
  
  ggbeeswarm::geom_quasirandom(groupOnX=FALSE, size = 4) +
  #geom_point() +
  facet_wrap(~set) +
  scale_color_manual(values = pals::brewer.set1(n = 10) %>% unname()) +
  scale_fill_manual(values = c('white','gray')) +
  cowplot::theme_cowplot() +
  ylab('Number of\nHVGs')
