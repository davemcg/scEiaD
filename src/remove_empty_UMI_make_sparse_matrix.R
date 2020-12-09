#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
save(args, file = 'testing/nu_reumi_args.Rdata')

#base_dir = args[5]
#SRS = args[1]
#REF = args[2]
outdir <- args[1]
mito_genelist <-scan(args[2], what = character(), sep = '\n')
srs_directories <- args[-(1:2)]
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

library(Seurat)
library(BUSpaRse)
library(tidyverse)
library(Matrix)
library(DropletUtils)
library(readr)
library(zeallot)

# input data from project
## its embarssing i didnt think of this first.


read_bus <- function(indir){
  sample_id <- str_extract(indir, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')
  l <- read_velocity_output(spliced_dir = indir,
                       spliced_name = "spliced",
                       unspliced_dir = indir,
                       unspliced_name = "unspliced") 
  colnames(l$spliced) <- paste(colnames(l$spliced), sample_id, sep = ':')
  if( length(colnames(l$spliced)) != length(colnames(l$spliced))  ){
    message(' more cells than barcodes')
    }
  
  colnames(l$unspliced) <- paste(colnames(l$unspliced), sample_id, sep = ':')
  return(l)
}

calc_spl_ratio <- function(x, cells){
  spl_sum <- sum(x$spliced[,cells])
  uns_sum <- sum(x$unspliced[,cells])
  return(uns_sum / (uns_sum + spl_sum))
  
}

remove_empty_droplets <- function(x, srs, mito_genelist){

  tot_count_spliced <- colSums(x$spliced)
  tot_count_unspliced <- colSums(x$unspliced)
  
  prefilter_common <- intersect(colnames(x$spliced), colnames(x$unspliced))
  
  bc_spliced <- tryCatch( expr = barcodeRanks(x$spliced), error=function(cond) return(NULL) )
  if(is.null(bc_spliced)){
    message('sample failed to meet umi threshold')
    df <- tibble(
      sample=srs,
      ncell_pf_spliced = -1,
      ncell_pf_unspliced = -1,
      pref_common  = -1,
      nbc_pass_spliced =  -1,
      nbc_pass_unspliced = -1,
      nbc_common = -1,
      nbc_union = -1,
      common_spl_ratio = -1,
      union_spl_ratio = -1
    )
    return(list(spliced = NULL,
                unspliced = NULL,
                stats = df))
  }
  bc_unspliced <- tryCatch( expr = barcodeRanks(x$unspliced), error=function(cond) return(NULL))
  if(is.null(bc_unspliced)){
    # keep only the spliced cells that meet filtering, but keep the other spliced cells as well,
    # to keep data the same shape
    bc_spliced_pass <- colnames(x$spliced)[tot_count_spliced > metadata(bc_spliced)$inflection]
    bc_unspliced_pass <- bc_spliced_pass
    df <- tibble(
      ncell_pf_spliced = ncol(x$spliced),
      ncell_pf_unspliced = ncol(x$unspliced),
      pref_common  = length(prefilter_common),
      nbc_pass_spliced = length(bc_spliced_pass),
      nbc_pass_unspliced = 0,
      nbc_common =-1,
      nbc_union = -1,
      common_spl_ratio = -1,
      union_spl_ratio = -1,
      
    )
  } else{
    bc_spliced_pass <- colnames(x$spliced)[tot_count_spliced > metadata(bc_spliced)$inflection]
    bc_unspliced_pass <- colnames(x$unspliced)[tot_count_unspliced > metadata(bc_unspliced)$inflection]
    
    
    common <- intersect(bc_spliced_pass, bc_unspliced_pass)
    bc_union <- union(bc_spliced_pass, bc_unspliced_pass) %>% intersect(prefilter_common)
    df <- tibble(
      ncell_pf_spliced = ncol(x$spliced),
      ncell_pf_unspliced = ncol(x$unspliced),
      pref_common  = length(prefilter_common),
      nbc_pass_spliced = length(bc_spliced_pass),
      nbc_pass_unspliced = length(bc_unspliced_pass),
      nbc_common = common %>% length,
      nbc_union = length(bc_union),
      common_spl_ratio = calc_spl_ratio(x, common),
      union_spl_ratio = calc_spl_ratio(x, bc_union),
      
    )
    
  }
  #keep cells that meet criteria in splcied data 
  common <- bc_spliced_pass %>% intersect(prefilter_common)
  spliced = x$spliced[,common] 
  unspliced = x$unspliced[,common]
  
  ## quality control: remove high mito cells, and remove high count(doublet) cells 
  seu <- CreateSeuratObject(spliced)
  cells_above_min_umi <-  seu$nFeature_RNA > 200
  cells_below_max_umi <- seu$nFeature_RNA < 3000 
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genelist)
  cells_below_max_mito_pt <-  seu$percent.mt < 10
  keep_cells <- cells_below_max_mito_pt & cells_above_min_umi & cells_below_max_umi
  df <- df %>% mutate(ncells_pre_qc = ncol(spliced), 
                      ncells_failed_min_umi = sum(!cells_above_min_umi), 
                      ncells_failed_max_umi = sum(!cells_below_max_umi), 
                      ncells_failed_mito = sum(!cells_below_max_mito_pt), 
                      ncells_total_pass_qc = sum(keep_cells))
  return(list(spliced = spliced[,keep_cells],
              unspliced = unspliced[,keep_cells],
              stats = df))

}


study_counts_list <- lapply(srs_directories,read_bus) 

srs_names <- str_extract(srs_directories, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')
names(study_counts_list) <- srs_names

filtered_counts <- lapply(seq_along(study_counts_list), function(i) remove_empty_droplets(study_counts_list[[i]], 
                                                                                         names(study_counts_list)[i],
                                                                                         mito_genelist))


spliced <- lapply(filtered_counts, function(x) x[['spliced']]) %>% purrr::reduce(RowMergeSparseMatrices)
unspliced <- lapply(filtered_counts, function(x) x[['unspliced']]) %>% purrr::reduce(RowMergeSparseMatrices)
stats <-  lapply(filtered_counts, function(x) x[['stats']]) %>% bind_rows

write_tsv(stats, path = stats_file)

# save pared down counts
save(spliced, file = spliced_matrix_file)
save(unspliced, file = unspliced_matrix_file)
message('finished successfully')

