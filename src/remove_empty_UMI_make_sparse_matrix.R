#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
base_dir <- '/data/mcgaugheyd/projects/nei/mcgaughey/massive_integrated_eye_scRNA/'

SRS = args[1]
REF = args[2]
matrix_file <- paste0(base_dir, args[3])
# tsne_plot <- args[3]
stats_file <- paste0(base_dir, args[4])

library(Seurat)
library(BUSpaRse)
library(Matrix)
library(DropletUtils)
library(readr)

# input data from project

raw_matrix <- BUSpaRse::read_count_output(paste0(base_dir, 'quant/', SRS, '/', REF, '/genecount'),'gene', FALSE)
dim(raw_matrix)
tot_counts <- Matrix::colSums(raw_matrix)


bc_rank <- try({ barcodeRanks(raw_matrix) })
if (class(bc_rank) == 'try-error') {
  bc_rank <- barcodeRanks(raw_matrix, lower = 50)
}
# qplot(bc_rank$total, bc_rank$rank, geom = "line") +
#   geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
#   geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
#   annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
#   label = c("knee", "inflection"), color = c("blue", "green")) +
#   scale_x_log10() +
#   scale_y_log10() +
#   labs(y = "Barcode rank", x = "Total UMI count")

res_matrix <- raw_matrix[, tot_counts > metadata(bc_rank)$inflection]
# dim(res_matrix)
# 
# seu <- CreateSeuratObject(res_matrix, min.cells = 3) %>%
#   NormalizeData(verbose = FALSE) %>%
#   ScaleData(verbose = FALSE) %>%
#   FindVariableFeatures(verbose = FALSE)
# 
# 
# seu <- RunPCA(seu, verbose = FALSE, npcs = 30)
# ElbowPlot(seu, ndims = 30)
# DimPlot(seu, reduction = "pca", pt.size = 0.5)
# 
# seu <- RunTSNE(seu, dims = 1:20, check_duplicates = FALSE)
# DimPlot(seu, reduction = "tsne", pt.size = 0.5)

# write out pre/post UMI counts
stats <- data.frame('Gene_Number' = c(dim(raw_matrix)[1], dim(res_matrix)[1]), 
                    'UMI_Count' = c(dim(raw_matrix)[2], dim(res_matrix)[2]),
                    'State' = c('Raw', 'Processed'),
                    'SRS' = c(SRS,SRS))
write_tsv(stats, path = stats_file)

# save pared down counts
save(res_matrix, file = matrix_file)
