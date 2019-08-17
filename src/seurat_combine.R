args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
plan(strategy = "multicore", workers = 6)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 2400000 * 1024^2)

# anchors
load(args[2])
seurat_merged <- IntegrateData(anchorset = anchors, normalization.method = 'SCT', verbose = TRUE)

# PCA, UMAP
seurat_merged <- RunPCA(seurat_merged, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:75, min.dist = 0.75)

# jackstraw calc
seurat_merged <- JackStraw(seurat_merged, num.replicate = 100, dims = 100)
#seurat_merged <- ScoreJackStraw(seurat_merged, dims = 1:100)

# clustering
seurat_merged <- FindNeighbors(seurat_merged,  reduction = "pca", dims = 1:75, nn.eps = 0.5)
seurat_merged <- FindClusters(seurat_merged, resolution = c(0.5), n.start = 50)

# find all markers 
markers <- FindAllMarkers(seurat_merged, min.pct = 0.5, max.cells.per.ident = 2000)

#pdf('merged.pdf')
#print(DimPlot( seurat_merged, group.by = 'orig.ident'))
#dev.off()

save(seurat_merged, markers, file = args[1], compress = FALSE)
