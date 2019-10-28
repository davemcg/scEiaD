load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__scran__downsample__batch__fastMNN__dims50.umap.Rdata')
fastmnn <- umap # best so far
load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__scran__downsample__batch__none__dims25.umap.Rdata')
none <- umap # poor integration
load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__standard__downsample__batch__harmony__dims25.umap.Rdata')
harmony <- umap # integrates well at study level, but scrambles cell identity




umap <- bind_rows(fastmnn, none, harmony)
# fix cell type coloring
cell_types <- umap %>% mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets', 'Fibroblasts', 'Red Blood Cells'),
         !grepl('RPE|Vascul', CellType)) %>% 
  pull(CellType) %>% unique() %>% sort()
type_val <- setNames(pals::alphabet(n = cell_types %>% length()), cell_types)
type_col <- scale_colour_manual(values = type_val, name = 'Cell Type')
type_fill <- scale_fill_manual(values = type_val, name = 'Cell Type')


pdf(height = 6, width = 10, 'figure1_integration_performance.pdf')
umap %>% filter(Age > 10) %>% 
  mutate(CellType_predict = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType_predict),
         Method = factor(Method, levels = c('harmony', 'none', 'fastMNN'))) %>% 
  filter(!is.na(CellType_predict), !is.na(organism), 
         CellType_predict %in% c('Amacrine Cells', 'Bipolar Cells', 'Muller Glia', 'Rods', 'Cones', 'Retinal Ganglion Cells')) %>% 
  ggplot(aes(x=UMAP_1, y = UMAP_2, colour = CellType_predict)) + 
  geom_point(size = 0.3, alpha = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + 
  facet_grid(vars(Method), vars(organism)) + 
  theme_cowplot() + 
  type_col + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + xlab('UMAP 1') + ylab('UMAP 2')
dev.off()
