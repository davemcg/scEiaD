#pseudotime assessment
library(tidyverse)
library(cowplot)

load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/monocle_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.monocle.Rdata')

load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.umap.Rdata')



umap <- left_join(umap, 
                  cds_retina@principal_graph_aux$UMAP$pseudotime %>% 
                    enframe(name = 'Barcode', value = 'pseudotime'))


# correlation
cor(umap %>% filter(Age < 500) %>% 
      filter(!is.na(pseudotime),
             pseudotime != Inf,
             !is.na(Age)) %>% pull(Age), 
    umap %>% filter(Age < 500) %>% 
      filter(!is.na(pseudotime),
             pseudotime != Inf,
             !is.na(Age)) %>% 
      pull(pseudotime), method = 'pearson')

# correlation
cor(umap %>% filter(Age < 20) %>% 
      filter(!is.na(pseudotime),
             pseudotime != Inf,
             !is.na(Age)) %>% pull(Age), 
    umap %>% filter(Age < 20) %>% 
      filter(!is.na(pseudotime),
             pseudotime != Inf,
             !is.na(Age)) %>% 
      pull(pseudotime), method = 'pearson')

umap %>% 
  filter(!is.na(pseudotime),
         pseudotime != Inf,
         !is.na(Age)) %>% 
  ggplot(aes(x=Age, y = pseudotime)) + 
  geom_boxplot(aes(group = Age)) +
  geom_smooth(method = 'lm') + 
  theme_cowplot()

umap %>% filter(Age < 20) %>% 
  filter(!is.na(pseudotime),
         pseudotime != Inf,
         !is.na(Age)) %>% 
  ggplot(aes(x=Age, y = pseudotime)) + 
  geom_boxplot(aes(group = Age)) +
  geom_smooth(method = 'lm') + 
  theme_cowplot()
