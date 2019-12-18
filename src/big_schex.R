library(schex)
library(tidyverse)
library(SingleCellExperiment)

cell_types <- umap %>% 
  pull(CellType_predict) %>% unique() %>% sort()
type_val <- setNames(pals::alphabet(n = cell_types %>% length()), cell_types)
type_col <- scale_colour_manual(values = type_val)
type_fill <- scale_fill_manual(values = type_val)


str(umap)

sce <- SingleCellExperiment(assay = list(counts = umap %>% select(RHO:NTNG1) %>% as.matrix() %>% t()))
reducedDims(sce) <- list(UMAP = umap %>% select(UMAP_1, UMAP_2) %>% as.matrix())
sce$CellType <- umap$CellType_predict
sce <- make_hexbin(sce, nbins = 150)
# hexbin celltype
plot_hexbin_meta(sce,  col = 'CellType', action = 'majority') + type_fill + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('Cell Type')

# rho
plot_hexbin_gene(sce, gene = 'RHO',type = 'counts', action = 'mean')

# opn1sw
plot_hexbin_gene(sce, gene = 'OPN1SW',type = 'counts', action = 'mean')




