# conos
library(Seurat)
library(SeuratWrappers)
library(conos)
library(tidyverse)
load('quant/Mus_musculus/scTransformCCA_merged_01.seuratV3.Rdata')
temp <- seurat_merged
DefaultAssay(temp) <- 'RNA'
panel <- SplitObject(temp, split.by = "orig.ident")
for (i in 1:length(panel)) {
  panel[[i]] <- NormalizeData(panel[[i]]) %>% FindVariableFeatures() %>% ScaleData() %>% 
    RunPCA(verbose = FALSE)
}
con <- Conos$new(panel, n.cores = 2)
con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 30, n.odgenes = 2000, matching.method = "mNN", 
                       metric = "angular", score.component.variance = TRUE, verbose = TRUE)
con$findCommunities()
con$embedGraph()
conos <- as.Seurat(con)

DimPlot(conos, reduction = "largeVis", group.by = c("orig.ident", "ident"), ncol = 3)