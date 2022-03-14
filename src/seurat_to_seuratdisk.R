library(Seurat)
library(SeuratDisk)

load('site/scEiaD_all_seurat_v3.Rdata')

SaveH5Seurat(scEiaD, filename='site/scEiaD_all_seurat_v3.h5Seurat', overwrite = TRUE)

