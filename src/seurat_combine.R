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
if (class(anchors) != 'list'){
  seurat_merged <- IntegrateData(anchorset = anchors, normalization.method = 'SCT', verbose = TRUE)
} else {
  seurat_merged <- list()
  for (i in names(anchors)){
    seurat_merged[[i]] <- IntegrateData(anchorset = anchors[[i]], normalization.method = 'SCT', verbose = TRUE)
  }
}

# PCA, UMAP
seurat_merged <- RunPCA(seurat_merged, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:75, min.dist = 0.2)

# jackstraw calc
seurat_merged <- JackStraw(seurat_merged, num.replicate = 100, dims = 100)
#seurat_merged <- ScoreJackStraw(seurat_merged, dims = 1:100)

# clustering
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:75, nn.eps = 0.5)
seurat_merged <- FindClusters(seurat_merged, 
                              resolution = c(0.1,0.3,0.6,0.8,1,2,3,4,5),
                              save.SNN = TRUE,
                              do.sparse = TRUE,
                              algorithm = 2,
                              random.seed = 23)

# clustree
library(clustree)
pdf('clustree_seurat.pdf', width = 30, height = 20)
clustree(seurat_merged)
dev.off()

seurat_merged <- FindClusters(seurat_merged, 
                              resolution = 1,
                              save.SNN = TRUE,
                              do.sparse = TRUE,
                              algorithm = 2,
                              random.seed = 23)

# find all markers 
markers <- FindAllMarkers(seurat_merged, min.pct = 0.5, max.cells.per.ident = 2000)

#pdf('merged.pdf')
#print(DimPlot( seurat_merged, group.by = 'orig.ident'))
#dev.off()

save(seurat_merged, markers, file = args[1], compress = FALSE)
