library(tidyverse)
library(singleCellHaystack)
library(Seurat)
library(Matrix)
library(tictoc)
args = commandArgs(trailingOnly=TRUE)
load(args[1])
num_cell <- as.integer(args[2])
set.seed(1253)
#if (is.na(num_cell)){num_cell = ncol(integrated_obj@assays$RNA@counts)}
#cells <- sample(1:ncol(integrated_obj@assays$RNA@counts), num_cell)
#detect <- ( integrated_obj@assays$RNA@counts[,cells] > 0)

if (!is.na(num_cell)){
	cells <- sample(1:ncol(integrated_obj@assays$RNA@counts), num_cell)
	cts <- integrated_obj@assays$RNA@counts[,cells]
	detect <- cts > 0
	gd <- colSums(detect)
	scvi <- Embeddings(integrated_obj, reduction = 'scVI')[cells,]
} else {
	detect <- ( integrated_obj@assays$RNA@counts > 0)
	scvi <- Embeddings(integrated_obj, reduction = 'scVI')#[cells,]
	cts <- ( integrated_obj@assays$RNA@counts) #[,cells])
	gd <- colSums(detect)
}
rm(integrated_obj)
tic()
scH <- haystack(x = scvi , detection = detect, use.advanced.sampling = gd)
toc()

save(scH, file = args[3])
