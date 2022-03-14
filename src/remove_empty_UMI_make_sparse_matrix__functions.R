# R

read_bus <- function(indir){
  sample_id <- str_extract(indir, glue('({patterns})\\d+') )
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

remove_empty_droplets <- function(x, srs, mito_genelist, v3, cells_above_min_umi_val = 200){
  print(paste0('cells_above_min_umi_val is: ', cells_above_min_umi_val))
  if (!is.na(v3)){
		mito_cutoff <- 20
	} else { mito_cutoff <- 10 
  }
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
      sample=srs,
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
      sample=srs,
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
  if(all(grepl('\\.\\d+$', sample(rownames(spliced), 10))) &
     all(!grepl('\\.\\d+$', mito_genelist, 10))) {
    #if rownames have version, but  mito genes do not, remove  versions before making seurat object
    s <- spliced# however, do not change the output rownames. the versions will be removed in the
    # next step
    rownames(s) <- str_remove_all(rownames(spliced), '\\.\\d+$')
    seu <- CreateSeuratObject(s)
  }else {
    seu <- CreateSeuratObject(spliced)
  }

  ## quality control: remove high mito cells, and remove high count(doublet) cells

  cells_above_min_umi <-  seu$nFeature_RNA > cells_above_min_umi_val
  cells_above_min_umi2 <-  seu$nFeature_RNA > 600
  cells_below_max_umi <- seu$nFeature_RNA < 100000

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genelist)
  pct_mt_df = tibble(srs=srs,barcode = colnames(seu), `percent.mt` = seu$percent.mt, nFeature_RNA = seu$nFeature_RNA)
  cells_below_max_mito_pt <-  seu$percent.mt < 20 
  keep_cells <- cells_below_max_mito_pt & cells_above_min_umi & cells_below_max_umi
  df <- df %>% mutate(ncells_pre_qc = ncol(spliced),
                      ncells_failed_min_umi = sum(!cells_above_min_umi),
                      ncells_failed_min_umi600 = sum(!cells_above_min_umi2),
                      ncells_failed_max_umi = sum(!cells_below_max_umi),
                      ncells_failed_mito = sum(!cells_below_max_mito_pt),
                      ncells_total_pass_qc = sum(keep_cells))
  if(sum(keep_cells)>1){
    spliced_out = spliced[,keep_cells]
    unspliced_out = unspliced[,keep_cells]
    stats_out = df
    pct_mt_df_out = pct_mt_df
  } else{# rare edge case where sample has only one cell

    spliced_out = Matrix(spliced[,keep_cells], ncol=1)
    unspliced_out = Matrix(unspliced[,keep_cells], ncol=1)
    colnames(spliced_out) <- colnames(unspliced_out) <- colnames(spliced)[keep_cells]
    rownames(spliced_out) <- rownames(spliced)
    rownames(unspliced_out) <- rownames(unspliced)
    stats_out = df
    pct_mt_df_out = pct_mt_df
  }
  return(list(spliced = spliced_out,
              unspliced = unspliced_out,
              stats = stats_out,
              pct_mt_df = pct_mt_df_out))

}

run_decontX <- function(filtered_count_matrix, unfiltered_count_matrix){
    sample = str_extract(colnames(filtered_count_matrix$spliced)[1], ':\\w+') %>% gsub(':','',.)
    print(sample)
    decontX_spliced <- try({decontX(x = filtered_count_matrix$spliced, background = unfiltered_count_matrix$spliced)})
    decontX_unspliced <- try({decontX(x = filtered_count_matrix$unspliced, background = unfiltered_count_matrix$unspliced)})
    if (class(decontX_spliced) == 'try-error'){
        decontX_spliced <- list()
        decontX_spliced$decontXcounts <- filtered_count_matrix$spliced
        print('decontX failed')
    } else {
        # write contmination stats
        # probably doublets?
        print('worked')
        contamination_spliced = decontX_spliced$contamination %>% enframe()
        contamination_spliced$barcode <- colnames(filtered_count_matrix$spliced)
        write_tsv(contamination_spliced, file = paste(trimws(outdir), paste0(trimws(sample), '_matrix_decontXcontamination_spliced.tsv.gz'), sep = '/'))
    }

    if (class(decontX_unspliced) == 'try-error'){
        decontX_unspliced <- list()
        decontX_unspliced$decontXcounts <- filtered_count_matrix$unspliced
        print('decontX failed')
    } else {
        contamination_unspliced = decontX_unspliced$contamination %>% enframe()
        contamination_unspliced$barcode <- colnames(filtered_count_matrix$unspliced)
        write_tsv(contamination_unspliced, file = paste(trimws(outdir), paste0((sample), '_matrix_decontXcontamination_unspliced.tsv.gz'), sep = '/'))
    }


    return(list(spliced = decontX_spliced$decontXcounts,
                unspliced = decontX_unspliced$decontXcounts
                ))
}


