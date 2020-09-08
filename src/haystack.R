library(tidyverse)
library(singleCellHaystack)
library(Seurat)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
load(args[1])
num_cell <- as.integer(args[2])
set.seed(1253)
cells <- sample(1:ncol(integrated_obj@assays$RNA@counts), num_cell)
detect <- as.matrix( integrated_obj@assays$RNA@counts[,cells] > 0)
scvi <- Embeddings(integrated_obj, reduction = 'scVI')[cells,]
cts <- as.matrix( integrated_obj@assays$RNA@counts[,cells])
gd <- apply(detect, 2, sum)
rm(integrated_obj)
scH <- haystack(x = scvi , detection = detect, use.advanced.sampling = gd)

save(scH, file = args[3])
