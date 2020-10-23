#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
system('mkdir -p testing')

base_dir = args[5]
SRS = args[1]
REF = args[2]
matrix_file_dir <- args[3]
stats_file <- args[4]

####
# base_dir = '/data/swamyvs/scEiaD/'
# SRS = 'SRS6424747'
# REF = 'hs-homo_sapiens'
# matrix_file_dir <- '/data/OGVFB_BG/new_quant_sciad/quant/SRS6424747/10xv2/hs-homo_sapiens/genecount/'
# stats_file <- args[4]
####


spliced_matrix_file <- paste(matrix_file_dir, 'matrix.Rdata', sep = '/')
unspliced_matrix_file = paste(matrix_file_dir, 'unspliced_matrix.Rdata', sep = '/')

library(Seurat)
library(BUSpaRse)
library(tidyverse)
library(Matrix)
library(DropletUtils)
library(readr)
library(zeallot)

# input data from project

c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = matrix_file_dir,
                                                spliced_name = "spliced",
                                                unspliced_dir = matrix_file_dir,
                                                unspliced_name = "unspliced")
tot_count <- Matrix::colSums(spliced)


bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)
bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]
# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 0]
sf <- spliced[genes_use, bcs_use]
uf <- unspliced[genes_use, bcs_use]


# write out pre/post UMI counts
stats <- tibble(pre_spliced_umi = sum(spliced) , 
                pre_unspliced_umi = sum(unspliced),
                post_spliced_umi = sum(sf),
                post_unspliced_umi = sum(uf),
                pre_spliced_gene_number = nrow(spliced) , 
                pre_unspliced_gene_number = nrow(unspliced),
                post_spliced_gene_number = nrow(sf),
                post_unspliced_gene_number = nrow(uf),
                pre_pt_uns  = sum(unspliced) /(sum(unspliced) +sum(spliced) ),
                post_pt_uns = sum(uf) /(sum(uf) +sum(sf) ),
                pt_uns_diff =   post_pt_uns - pre_pt_uns
)
print(stats)
write_tsv(stats, path = stats_file)

# save pared down counts
save(sf, file = spliced_matrix_file)
save(uf, file = unspliced_matrix_file)
message('finished successfully')

