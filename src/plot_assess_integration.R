library(tidyverse)
library(Polychrome)
args <- commandArgs(trailingOnly = TRUE)

load(args[1])

# color by study, facet by age
pdf(args[2], height = 10, width = 12)
umap %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = study_accession)) + 
  geom_point(size = 0.3, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~Age) + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()

# color by study_accession, facet by batch
pdf(args[3], height = 10, width = 12)
umap %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = study_accession)) + 
  geom_point(size = 0.3, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~batch) + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()

# color by known labelling paper, facet by celltype
pdf(args[4], height = 10, width = 12)
umap %>% mutate(Labelling = case_when(is.na(Paper) ~ 'None',
                                      TRUE ~ Paper)) %>% 
  filter(!CellType %in% c('Doublets', 'Red Blood Cells')) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = Labelling)) + 
  geom_point(size = 0.3, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~CellType) + 
  theme_minimal() + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()

# color by celltype, facet by seurat cluster
pdf(args[5], height = 10, width = 12)
umap %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = CellType)) + 
  geom_point(size = 0.6, alpha = 0.1) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~seurat_clusters) + 
  scale_color_manual(values = unname(alphabet.colors())) +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()
