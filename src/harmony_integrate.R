rdata_files <- c('quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

setwd('/data/mcgaugheyd/projects/nei/mcgaughey/massive_integrated_eye_scRNA')
library(Seurat)
library(harmony)
library(tidyverse)

temp <- seurat_merged
DefaultAssay(temp) <- 'RNA'
harmony <- NormalizeData(temp) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
# crucial parameter is iterations. More iterations bring closer together, but can overfit?
harmony <- RunHarmony(harmony, group.by.vars = "orig.ident",
                      max.iter.harmony = 20, 
                      epsilon.harmony = -Inf)
harmony <- RunUMAP(harmony, 
                   reduction = "harmony", 
                   dims = 1:30)
harmony <- RunUMAP(harmony, 
                   n.components=3, 
                   reduction.name = 'umap3D', 
                   reduction = "harmony", 
                   dims = 1:30)
harmony <- FindNeighbors(harmony, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = c(0.1,0.3,0.6,0.8,1,2),
               save.SNN = TRUE,
               do.sparse = TRUE,
               algorithm = 2,
               random.seed = 23)
# run clustree to pick resolution
# clustree(harmony)
# 0.8 looks good (which is the default....)
harmony <- FindClusters(harmony, resolution = 0.8)

harmony_markers <- FindAllMarkers(harmony, 
                                  logfc.threshold=0.5, 
                                  max.cells.per.ident=3000, 
                                  min.pct=0.25)

save(harmony, harmony_markers, file = 'harmony_mouse.Rdata', compress = FALSE)