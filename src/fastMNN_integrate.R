#fastMNN
library(Seurat)
library(SeuratWrappers)
library(clustree)
load('quant/Mus_musculus/scTransformCCA_merged_01.seuratV3.Rdata')

temp <- seurat_merged
# reset defaultAssay to 'RNA' from 'integrated'
DefaultAssay(temp) <- 'RNA'
temp <- FindVariableFeatures(temp)
temp <- RunFastMNN(object.list = SplitObject(temp, split.by = "orig.ident"))
temp <- RunUMAP(temp, reduction = "mnn", dims = 1:30)
temp <- FindNeighbors(temp, reduction = "mnn", dims = 1:30)
temp <- FindClusters(temp, resolution = c(0.1,0.3,0.6,0.8,1,2,3,4,5),
                     save.SNN = TRUE,
                     do.sparse = TRUE,
                     algorithm = 2,
                     random.seed = 23)
mnn <- temp
# run clustree to pick resolution
pdf('mnn_clustree.pdf', width = 10, height = 15)
clustree(mnn)
dev.off()
# 2 looks good?
mnn <- FindClusters(mnn, resolution = 2)

mnn_markers <- FindAllMarkers(mnn, 
                                  logfc.threshold=0.5, 
                                  max.cells.per.ident=3000, 
                                  min.pct=0.25)

save(mnn, mnn_markers, file = args[2], compress = FALSE)