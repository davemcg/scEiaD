library(tidyverse)
library(Polychrome)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])

# specific plot with Bipolar Cells, Muller Glia, and Rods for >10 day old samples
# these are three cell types with multiple studies (matched roughly by time)
# would expect them to largely overlap (within each cell type)
# in other words sholud overlap wihtin each plot, but be in distinct space (x - y)
# for each cell type plot facet
pdf(args[2], height = 6, width = 8)
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
dev.odf()