library(tidyverse)
library(Polychrome)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])
colnames(umap) <- gsub('^rna_', '', colnames(umap))
# specific plot with Bipolar Cells, Muller Glia, and Rods for >10 day old samples
# these are three cell types with multiple studies (matched roughly by time)
# would expect them to largely overlap (within each cell type)
# in other words sholud overlap wihtin each plot, but be in distinct space (x - y)
# for each cell type plot facet
pdf(args[2], height = 4, width = 10)
umap %>% filter(Age > 10) %>% mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), Labelling = case_when(is.na(Paper) ~ 'None',
                                                                                                                            TRUE ~ Paper)) %>% 
  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods')) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = Labelling)) + 
  geom_point(size = 0.3, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~CellType) + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()

pdf(args[3], height = 4, width = 10)
p1 <- umap %>% filter(Age > 10) %>% 
  mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
         Labelling = case_when(is.na(Paper) ~ 'None',
                               TRUE ~ Paper)) %>% 
  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods')) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = CellType)) + 
  geom_point(size = 0.3, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p2 <- bind_rows(umap %>% filter(Age > 10) %>% mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
                                                     Labelling = case_when(is.na(Paper) ~ 'None',
                                                                           TRUE ~ Paper)) %>% 
                  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods')) %>% 
                  filter(CABP5 > 1) %>% mutate(Expression = 'Cabp5 (Bipolar Cell)'),
                
                umap %>% filter(Age > 10) %>% mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
                                                     Labelling = case_when(is.na(Paper) ~ 'None',
                                                                           TRUE ~ Paper)) %>% 
                  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods')) %>% 
                  filter(RHO > 2) %>% mutate(Expression = 'Rho (Rod Photoreceptor)'),
                
                umap %>% filter(Age > 10) %>% mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
                                                     Labelling = case_when(is.na(Paper) ~ 'None',
                                                                           TRUE ~ Paper)) %>% 
                  filter(CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods')) %>% 
                  filter(AQP4 > 1) %>% mutate(Expression = 'Aqp4 (Muller Glia)')) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = Expression)) + 
  geom_point(size = 0.3, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  scale_color_brewer(palette = 'Set2') + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
plot_grid(p1, p2)
dev.off()


