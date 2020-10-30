#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
save(args, file = 'testing/nu_reumi_args.Rdata')

#base_dir = args[5]
#SRS = args[1]
#REF = args[2]
outdir <- args[1]

srs_directories <- args[-1]
####
# base_dir = '/data/swamyvs/scEiaD/'
# SRS = 'SRS6424747'
# REF = 'hs-homo_sapiens'
# matrix_file_dir <- '/data/OGVFB_BG/new_quant_sciad/quant/SRS6424747/10xv2/hs-homo_sapiens/genecount/'
# stats_file <- args[4]
####


spliced_matrix_file <- paste(outdir, 'matrix.Rdata', sep = '/')
unspliced_matrix_file = paste(outdir, 'unspliced_matrix.Rdata', sep = '/')
stats_file <- paste(outdir, 'stats.tsv', sep = '/')

library(Seurat)
library(BUSpaRse)
library(tidyverse)
library(Matrix)
library(DropletUtils)
library(readr)
library(zeallot)

# input data from project
## its embarssing i didnt think of this first.

## SRS5396948
read_bus <- function(indir){
  sample_id <- str_extract(indir, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')
  l <- read_velocity_output(spliced_dir = indir,
                       spliced_name = "spliced",
                       unspliced_dir = indir,
                       unspliced_name = "unspliced") 
  colnames(l$spliced) <- paste(colnames(l$spliced), sample_id, sep = ':')
  colnames(l$unspliced) <- paste(colnames(l$unspliced), sample_id, sep = ':')
  return(l)
}

study_counts_list <- lapply(srs_directories,read_bus) 
spliced <- lapply(study_counts_list, function(x) x[['spliced']]) %>% purrr::reduce(RowMergeSparseMatrices)
unspliced <- lapply(study_counts_list, function(x) x[['unspliced']]) %>% purrr::reduce(RowMergeSparseMatrices)


tot_count <- Matrix::colSums(spliced)


bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)
bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]
# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 0]
sf <- spliced[genes_use, bcs_use]
uf <- unspliced[genes_use, bcs_use[bcs_use %in% colnames(unspliced)]]


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

