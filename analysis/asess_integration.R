# celltype prediction accuracy
load('~/data/scEiaD_v3/merged_xgboost_2021-11-05.Rdata')
load('~/data/scEiaD_v3/merged_stats_2021-11-05.Rdata')
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


#perf <- perf %>% filter(grepl('scVI',method)) 
#cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 

#cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 
# 
# perf_well <- perf %>% unique() %>% filter(set == 'onlyWELL') %>% 
#   select(Score, Group, Value, set, dims:normalization) %>% 
#   filter(Score %in% c('LISI','Silhouette')) 

# lisi <- perf %>% 
#   filter(Score == 'LISI') %>%
#   pivot_wider(values_from = Value, names_from = c('Group')) %>% 
#   ggplot(aes(y=Batch, x=-Cluster)) + geom_point(aes(color=as.factor(dims), shape = as.factor(nf)), size = 5) +
#   cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
#   ggtitle('LISI') 
# 
# ari <- perf %>% 
#   filter(Score == 'ARI', Group == 'CellType-Cluster') %>%
#   ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=as.factor(dims), shape = as.factor(nf)), size = 5, width = 0.2) +
#   cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('ARI') +
#   ggtitle('ARI')
# 
# nmi <- perf %>% 
#   filter(Score == 'NMI', Group == 'CellType-Cluster', nf == 2000, dims %in% c(8,30), knn == 20) %>%
#   ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
#   cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('NMI') +
#   ggtitle('NMI') + theme(legend.position="none") +
#   scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))
# 
# silhouette <- perf %>% filter(Score == 'Silhouette') %>% 
#   pivot_wider(values_from = Value, names_from = c('Group')) %>% 
#   ggplot(aes(y=-Batch, x=Cluster, shape = as.factor(nf), color = as.factor(dims))) + geom_point(size = 5) +
#   cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
#   ggtitle("Silhouette")
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   lisi + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
accuracy %>% 
  filter(partition %in% c('universe')) %>% 
  group_by(norm, nf, dims, method, knn, epochs) %>% 
  summarise(bad_score_sum = sum(score < 0.8), xgboost_score = mean(score)) %>% 
  arrange(-xgboost_score)

purity <- 1; balance <- 3; mixing <- 1
partition = 'universe'
method = 'scvIprojection'
perf %>% 
  filter(partition %in% c('universe')) %>% 
  select(-clusterN, -clusterMedian, -subset) %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  left_join(accuracy %>% 
              filter(partition %in% c('universe')) %>% 
              group_by(norm, nf, dims, method, knn, epochs) %>% 
              summarise(bad_score_sum = sum(score < 0.8), xgboost_score = mean(score))) %>% 
  mutate(sumZScale =
           -scale(LISI_CellType)[,1] * purity +
           scale(Silhouette_CellType)[,1] * purity  +
           scale(Silhouette_Cluster)[,1] * purity  +
           -scale(LISI_Cluster)[,1] * purity  +
           
           scale(`xgboost_score`)[,1] * balance +
           
           scale(`NMI_CellType-Cluster`)[,1] * balance  +
           scale(`ARI_CellType-Cluster`)[,1]* balance  +
           
           scale(LISI_Batch)[,1] * mixing  %>%   # Z score
           -scale(Silhouette_Batch)[,1] * mixing) %>% 
  mutate(nf = as.factor(nf),
         `scVI latent dims` = as.factor(dims),
         method = gsub('scVIprojection', 'scVI-projection',method),
         method = gsub('scVI$', 'scVI-standard', method),
         method = factor(method, levels = c('scVI-standard','scVI-projection'))) %>% 
  ggplot(aes(x=sumZScale, y = nf, color = `scVI latent dims`, shape = as.factor(epochs))) + 
  ggbeeswarm::geom_quasirandom(groupOnX=FALSE, size = 5) +
  #geom_point() +
  facet_wrap(~set) +
  scale_color_manual(values = pals::brewer.set1(n = 10) %>% unname()) +
  scale_fill_manual(values = c('white','gray')) +
  cowplot::theme_cowplot() +
  ylab('Number of\nHVGs')
