#!/usr/bin/env Rscript/

args = commandArgs(trailingOnly=TRUE)
# save(args, file = 'testing/nu_reumi_args.Rdata')

outdir <- args[1]
mito_genelist <-scan(args[2], what = character(), sep = '\n')
# patterns <- args[3]
git_dir <- args[5]
srs_directories <- args[-(1:5)]
########################################################
# base_dir = '/data/swamyvs/scEiaD/'
# SRS = 'SRS6424747'
# REF = 'hs-homo_sapiens'
# matrix_file_dir <- '/data/OGVFB_BG/new_quant_sciad/quant/SRS6424747/10xv2/hs-homo_sapiens/genecount/'
# stats_file <- args[4]
################################################


spliced_matrix_file <- paste(outdir, 'matrix.Rdata', sep = '/')
unspliced_matrix_file = paste(outdir, 'unspliced_matrix.Rdata', sep = '/')
stats_file <- paste(outdir, 'stats.tsv', sep = '/')
pct_mt_file = paste(outdir, 'pct_mt.tsv', sep='/')
library(Seurat)
library(BUSpaRse)
library(tidyverse)
library(Matrix)
library(DropletUtils)
library(readr)
library(zeallot)
library(glue)
library(celda) #decontX
# input data from project
## its embarssing i didnt think of this first.
patterns <- scan(args[3], what = character(), sep='\n') %>% paste0(collapse = '|')
ambient <- args[4]

# load scripts
source(glue('{git_dir}/src/remove_empty_UMI_make_sparse_matrix__functions.R'))

# data import
study_counts_list <- lapply(srs_directories,read_bus) 

# grab srs names
srs_names <- str_extract(srs_directories, glue('({patterns})\\d+') )
names(study_counts_list) <- srs_names

# find 10xv3
v3 <- str_extract(srs_directories, '10xv3')
# run filtering
filtered_counts <- lapply(seq_along(study_counts_list), function(i) remove_empty_droplets(study_counts_list[[i]], 
                                                                                         names(study_counts_list)[i],
                                                                                         mito_genelist,
																						 v3[i]))
filtered_counts_orig <- filtered_counts
# if there is a sample with fewer than 50 cells, then just combine them all
low_n_count <- 0
for (i in seq_along(filtered_counts)){
	if ( (filtered_counts[[i]]$spliced %>% ncol() ) < 30){
		low_n_count <- low_n_count + 1 
	}
}
# run decontX ambient RNA removal
if (ambient == 'decontX'){
	library(parallel)
	decontX_counts <- mclapply(seq_along(study_counts_list), mc.cores = 4, function(i) run_decontX(filtered_counts[[i]],
																				study_counts_list[[i]]))
	# rerun remove_empty_droplets to remove cells dropping under 300 counts
	decontX_counts_stats <- lapply(seq_along(study_counts_list), function(i) remove_empty_droplets(decontX_counts[[i]],
                                                                                         names(study_counts_list)[i],
                                                                                         mito_genelist,
																						 v3[i]))
	# deal with edge case where no cells make it past the filtering (remove_empty_droplets) part 2
	for (i in seq_along(study_counts_list)){
		if (is.null(decontX_counts_stats[[i]]$spliced)){
			decontX_counts_stats[[i]]$spliced <- filtered_counts[[i]]$spliced
		}
		if (is.null(decontX_counts_stats[[i]]$unspliced)){
			decontX_counts_stats[[i]]$unspliced <- filtered_counts[[i]]$unspliced
		}
	}
	filtered_counts <- decontX_counts_stats
}
spliced <- lapply(filtered_counts, function(x) x[['spliced']]) %>% purrr::reduce(RowMergeSparseMatrices)
unspliced <- lapply(filtered_counts, function(x) x[['unspliced']]) %>% purrr::reduce(RowMergeSparseMatrices)
stats <-  lapply(filtered_counts_orig, function(x) x[['stats']]) %>% bind_rows
pct_mt = lapply(filtered_counts_orig, function(x) x[['pct_mt_df']]) %>% bind_rows
write_tsv(stats, path = stats_file)
write_tsv(pct_mt, path = pct_mt_file)
# save pared down counts
save(spliced, file = spliced_matrix_file)
save(unspliced, file = unspliced_matrix_file)
message('finished successfully')

