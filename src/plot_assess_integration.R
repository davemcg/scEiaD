library(tidyverse)
library(pals)
library(cowplot)
args <- commandArgs(trailingOnly = TRUE)

load(args[1])

# color by study, facet by age
pdf(args[2], height = 10, width = 12)
umap %>% 
  mutate(Time = integration_group,
         CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
         Labelling = case_when(Paper == 'Hufnagel 2020' ~ 'Homo sapiens (Stem) (Hufnagel 2020)',
                               TRUE ~ paste0(organism, ' (', study_accession, ')'))) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = Labelling)) + 
  geom_point(size = 0.1, alpha = 0.05) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~Age) + 
  theme_cowplot() + 
  scale_color_manual(values = as.vector(pals::alphabet())) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()

# color by study_accession, facet by batch
pdf(args[3], height = 10, width = 12)
umap %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = study_accession)) + 
  geom_point(size = 0.1, alpha = 0.05) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~batch) + 
  theme_cowplot() + 
  scale_color_manual(values = as.vector(pals::alphabet())) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
dev.off()

# color by known labelling paper, facet by celltype
pdf(args[4], height = 10, width = 12)
remove_low_n <- umap %>% mutate(Time = integration_group,
                                CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
                                Labelling = case_when(is.na(Paper) ~ 'None',
                                                      Paper == 'Hufnagel 2020' ~ 'Homo sapiens (Stem) (Hufnagel 2020)',
                                                      TRUE ~ paste0(organism, ' (', Paper, ')'))) %>% 
  filter(!is.na(CellType), !CellType %in% c('Doublet', 'Doublets', 
                                            'Red Blood Cells', 'Fibroblasts',
                                            'Astrocytes', 'Pericytes',
                                            'RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cells')) %>% group_by(CellType, Time) %>% 
  summarise(Count = n(), Barcodes = list(Barcode)) %>% filter(Count < 250) %>% pull(Barcodes) %>% unlist()
umap %>% mutate(Time = integration_group,
                CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType), 
                Labelling = case_when(is.na(Paper) ~ 'None',
                                      Paper == 'Hufnagel 2020' ~ 'Homo sapiens (Stem) (Hufnagel 2020)',
                                      TRUE ~ paste0(organism, ' (', Paper, ')'))) %>% 
  filter(!is.na(CellType), 
         !Barcode %in% remove_low_n, 
         !CellType %in% c('Doublet', 'Doublets', 
                                            'Red Blood Cells', 'Fibroblasts',
                                            'Astrocytes', 'Pericytes',
                                            'RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cells')) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = Labelling)) + 
  geom_point(size = 0.1, alpha = 0.05) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_wrap(~CellType + Time) + 
  theme_cowplot() + 
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
