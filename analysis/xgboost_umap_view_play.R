library(tidyverse)
library(cowplot)
library(scattermore)
library(glue)
load('~/data/scEiaD_v2/xgboost_predictions/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-15__epochs-5__dist-0.1__neighbors-50__knn-20__umapPredictions.Rdata')
ncells <- nrow(umap)

labels <- umap %>% group_by(CellType_predict) %>% 
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))
# pan view
umap %>% 
  mutate(CellType_predict = case_when(is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) %>% 
  ggplot() + 
  geom_scattermore(aes(x=UMAP_1, 
                       y = UMAP_2, 
                       colour = CellType_predict), pointsize = 0, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  ggrepel::geom_label_repel(data = labels, aes(x=UMAP_1, y=UMAP_2, label = CellType_predict ), alpha = 0.8, size = 3) +
  scale_color_manual(values = c(pals::alphabet(), pals::alphabet2(), pals::glasbey()) %>% unname()) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  ggtitle(glue('Number of Cells: {ncells}')) 


## cluster
umap %>% 
  ggplot() + 
  geom_scattermore(aes(x=UMAP_1, 
                       y = UMAP_2, 
                       colour = as.factor(cluster)), pointsize = 0, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() + 
  ggrepel::geom_label_repel(data = labels, aes(x=UMAP_1, y=UMAP_2, label = CellType_predict ), alpha = 0.8, size = 3) +
  scale_color_manual(values = c(pals::alphabet(), pals::alphabet2(), pals::glasbey()) %>% unname()) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  ggtitle(glue('Number of Cells: {ncells}')) 


# facet by cell type
umap %>% 
  filter(study_accession != 'SRP131661') %>% 
  mutate(CellType_predict = case_when(is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) %>% 
  ggplot() + 
  geom_scattermore(aes(x=UMAP_1, 
                       y = UMAP_2, 
                       colour = study_accession), pointsize = 5, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() +
  scale_color_manual(values = c(pals::alphabet2(), pals::alphabet(), pals::glasbey()) %>% unname()) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  ggtitle(glue('Number of Cells: {ncells}')) + facet_wrap(~CellType_predict) 

# RPE view
umap %>% filter(CellType_predict == 'RPE') %>% ggplot() + 
  geom_scattermore(aes(x=UMAP_1, 
                       y = UMAP_2), pointsize = 5, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=8, alpha = 1))) + 
  theme_cowplot() +
  scale_color_manual(values = c(pals::alphabet2(), pals::alphabet(), pals::glasbey()) %>% unname()) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  ggtitle(glue('Number of Cells: {ncells}')) + facet_wrap(~study_accession) 


# Quick Bharti RPE assess (should be mostly RPE/Fibro/Endo/Epi)
umap %>% group_by(CellType_predict, study_accession) %>% summarise(Count = n()) %>% filter(study_accession == 'Bharti_Nguyen_iRPE_2D_3D') %>% arrange(-Count)
# Cluster <-> CellType_predict
umap %>% 
  mutate(cluster = as.factor(cluster), CellType_predict = as.factor(CellType_predict)) %>% 
  group_by(cluster, CellType_predict, study_accession) %>% summarise(Count = n()) %>% mutate(Percent = (Count / sum(Count)) * 100) %>% filter(Percent > 5) %>% arrange(cluster, -Percent) %>% DT::datatable(
                                                                                                                                                                                                  filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                                                                                                                                                                                                  options = list(pageLength = 10,  searchHighlight = TRUE, dom = 'frtBip'))
