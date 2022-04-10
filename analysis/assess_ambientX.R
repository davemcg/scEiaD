library(tidyverse)
gene_name %>% filter(Name == 'RHO')
gene_name <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/gene_name_ids.tsv.gz')
load('~/data/scEiaD/ambientX/ambientX__n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-10__epochs-5__dist-0.1__neighbors-50__umap.seuratV3.Rdata')
ambientX_obj <- integrated_obj
load('~/data/scEiaD/ambientX/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-10__epochs-5__dist-0.1__neighbors-50__umap.seuratV3.Rdata')
a <- integrated_obj@assays$RNA@counts['ENSG00000185527',] %>% 
  enframe() %>% 
  rename(Barcode = name) %>% 
  left_join(umap) %>% 
  group_by(CellType_predict, organism) %>% 
  summarise(Mean = mean(log2(value+1) , na.rm = TRUE)) %>% 
  mutate(Core = case_when(CellType_predict == 'Rods' ~ 'True', TRUE ~ 'Contamination')) %>% 
  filter(!is.na(CellType_predict)) %>% 
  ggplot(aes(x=Core, y = Mean)) + 
  geom_boxplot(outlier.colour = NA) + 
  ggbeeswarm::geom_quasirandom(aes(color = organism), alpha = 0.5) + 
  cowplot::theme_cowplot() + 
  ggtitle('Raw Counts')

b <- ambientX_obj@assays$RNA@counts['ENSG00000185527',] %>% 
  enframe() %>% 
  rename(Barcode = name) %>% 
  left_join(umap) %>% 
  group_by(CellType_predict, organism) %>% 
  summarise(Mean = mean(log2(value+1) , na.rm = TRUE)) %>% 
  mutate(Core = case_when(CellType_predict == 'Rods' ~ 'True', TRUE ~ 'Contamination')) %>% 
  filter(!is.na(CellType_predict)) %>% 
  ggplot(aes(x=Core, y = Mean)) + 
  geom_boxplot(outlier.colour = NA) + 
  ggbeeswarm::geom_quasirandom(aes(color = organism), alpha = 0.5) + 
  cowplot::theme_cowplot() + 
  ggtitle('ambientX Counts')

plot_grid(a, b, ncol = 2)

# dumbbell
combined <- bind_rows(
  integrated_obj@assays$RNA@counts['ENSG00000163914',] %>% 
    enframe() %>% 
    rename(Barcode = name) %>% 
    left_join(umap) %>% 
    group_by(CellType_predict, organism) %>% 
    summarise(Mean = mean(log2(value+1) , na.rm = TRUE)) %>% 
    mutate(Counts = 'Raw'),
  ambientX_obj@assays$RNA@counts['ENSG00000163914',] %>% 
    enframe() %>% 
    rename(Barcode = name) %>% 
    left_join(umap) %>% 
    group_by(CellType_predict, organism) %>% 
    summarise(Mean = mean(log2(value+1) , na.rm = TRUE)) %>% 
    mutate(Counts = 'ambientX')
)

combined %>% 
  filter(!is.na(CellType_predict)) %>% 
  ggplot(aes(x=CellType_predict, y = Mean)) + 
  geom_line() + 
  coord_flip() + 
  facet_wrap(~organism) + 
  geom_point(aes(color = Counts)) + 
  cowplot::theme_cowplot() + 
  ylab("Mean Rhodopsin Expression \n(Rod Specific Marker)") +
  ggsci::scale_color_aaas()












combined2 <- bind_rows(
  integrated_obj@assays$RNA@counts['ENSG00000163914',] %>% 
    enframe() %>% 
    rename(Barcode = name) %>% 
    left_join(umap) %>% 
    #filter(TechType %in% c("10xv2",'10xv3','DropSeq')) %>% 
    group_by(CellType_predict, study_accession) %>% 
    summarise(Mean = mean(log2(value+1) , na.rm = TRUE)) %>% 
    mutate(Counts = 'Raw'),
  ambientX_obj@assays$RNA@counts['ENSG00000163914',] %>% 
    enframe() %>% 
    rename(Barcode = name) %>% 
    left_join(umap) %>% 
    #filter(TechType %in% c("10xv2",'10xv3','DropSeq')) %>% 
    group_by(CellType_predict, study_accession) %>% 
    summarise(Mean = mean(log2(value+1) , na.rm = TRUE)) %>% 
    mutate(Counts = 'ambientX')
)

combined2 %>% 
  filter(!is.na(CellType_predict)) %>% 
  ggplot(aes(x=CellType_predict, y = Mean)) + 
  geom_line() + 
  coord_flip() + 
  facet_wrap(~study_accession, ncol = 12) + 
  geom_point(aes(color = Counts)) + 
  cowplot::theme_cowplot() + 
  ylab("Mean Rhodopsin Expression \n(Rod Specific Marker)") +
  ggsci::scale_color_aaas()





combined2 %>% 
  filter(!is.na(CellType_predict), grepl('SRP21278|SRP21215|SRP269635', study_accession)) %>% 
  ggplot(aes(x=CellType_predict, y = Mean)) + 
  geom_line() + 
  coord_flip() + 
  facet_wrap(~study_accession, ncol = 12) + 
  geom_point(aes(color = Counts)) + 
  cowplot::theme_cowplot() + 
  ylab("Mean Rhodopsin Expression \n(Rod Specific Marker)") +
  ggsci::scale_color_aaas()

  