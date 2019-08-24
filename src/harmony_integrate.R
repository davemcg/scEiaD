args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(harmony)
library(tidyverse)

# load integrated Seurat object
load(args[1])
temp <- seurat_merged
# reset defaultAssay to 'RNA' from 'integrated'
DefaultAssay(temp) <- 'RNA'
# normalize for Harmony
harmony <- NormalizeData(temp) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
# crucial parameter is iterations. More iterations bring closer together, but can overfit?
harmony <- RunHarmony(harmony, group.by.vars = "orig.ident",
                      max.iter.harmony = 20, 
                      epsilon.harmony = -Inf)
harmony <- RunUMAP(harmony, 
                   reduction = "harmony", 
                   dims = 1:30)
# 3D umap may be helpful
harmony <- RunUMAP(harmony, 
                   n.components=3, 
                   reduction.name = 'umap3D', 
                   reduction = "harmony", reduction.key = 'umap3D',
                   dims = 1:30)
# run at many resolutions to pick best one (via clustree or cluster purity)
harmony <- FindNeighbors(harmony, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = c(0.1,0.3,0.6,0.8,1,2,3,4,5),
               save.SNN = TRUE,
               do.sparse = TRUE,
               algorithm = 2,
               random.seed = 23)

# run clustree to pick resolution
pdf('harmony_clustree.pdf', width = 10, height = 15)
clustree(harmony)
dev.off()
# 2 looks good?
harmony <- FindClusters(harmony, resolution = 2)

harmony_markers <- FindAllMarkers(harmony, 
                                  logfc.threshold=0.5, 
                                  max.cells.per.ident=3000, 
                                  min.pct=0.25)

save(harmony, harmony_markers, file = args[2], compress = FALSE)