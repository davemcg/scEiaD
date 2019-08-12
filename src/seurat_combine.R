args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)

# anchors
load(args[2])
seurat_merged <- IntegrateData(anchorset = anchors, normalization.method = 'SCT', verbose = TRUE)

# PCA, UMAP
seurat_merged <- RunPCA(seurat_merged, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:75, min.dist = 0.75)

# clustering
seurat_merged <- FindNeighbors(seurat_merged,  reduction = "pca", dims = 1:50, nn.eps = 0.5)
seurat_merged <- FindClusters(seurat_merged, resolution = c(0.5,1,2,3,5), n.start = 50)

#pdf('merged.pdf')
#print(DimPlot( seurat_merged, group.by = 'orig.ident'))
#dev.off()

save(seurat_merged, file = args[1])
