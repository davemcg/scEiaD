### 
# Use >= R/4.0
library(tictoc)
library(tidyverse)
library(scater)
library(Seurat)
library(edgeR)
library(BiocParallel)
multicoreParam <- MulticoreParam(workers = 12)
git_dir = Sys.getenv('SCIAD_GIT_DIR')
args <- commandArgs(trailingOnly = TRUE)

load(args[1]) #edgeR obj with fit, dispersions, and processed data
comp <- args[2]
partition = args[3] %>% as.numeric()
out <- args[4]
###############
# functions -------
###############

source( paste0(git_dir, '/src/pseudoBulk_functions.R') )
#####################

processed_data <- edgeR_obj$processed_data

if (grepl('A', comp)){
  ######################
  # celltype (pre-labelled/published) -------
  ####################
  
  #processed_data <- processing(sum_mat1)
  
  if (grepl('1', comp)){
  # celltype against remaining, controlling for organism ------
  CELLTYPE__res_againstAll <- pseudoBulk_testing(processed_data, 
                                                 organism_covariate=TRUE,
                                                 pairwise=FALSE,
												 edgeR_obj = edgeR_obj,
														partition = partition)
  save(CELLTYPE__res_againstAll, file = out)
  } else if (grepl('2', comp)){
  # celltype against each celltype (pairwise), controlling for organism ----------
  CELLTYPE__res_pairwise <- pseudoBulk_testing(processed_data, 
                                               organism_covariate=TRUE,
                                               pairwise=TRUE,
								               pieces = 500,
											   edgeR_obj = edgeR_obj,
														partition = partition)
  save(CELLTYPE__res_pairwise, file = out)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE
  CELLTYPE__res_organism_celltype <- try({pseudoBulk_testing(processed_data, 
                                                        organism_covariate=FALSE,
                                                        pairwise=TRUE, 
                                                        testing_against_internal_organism = TRUE,
                                                        testing_against = 'var_organism',
								                        edgeR_obj = edgeR_obj,
														partition = partition) })
  save(CELLTYPE__res_organism_celltype, file = out)
  }
} else if (grepl('B', comp)){
  ##############################
  # same, but with celltype predict (machine label all cells with label) ----------
  ##############################
  
  if (grepl('1', comp)){
  # CellType_predict against remaining, controlling for organism ------
  CELLTYPEPREDICT__res_againstAll <- pseudoBulk_testing(processed_data, 
                                                        organism_covariate=TRUE,
                                                        pairwise=FALSE,
												        edgeR_obj = edgeR_obj,
														partition = partition)
  save(CELLTYPEPREDICT__res_againstAll, file = out)
  } else if (grepl('2', comp)){
  # CellType_predict against each celltype (pairwise), controlling for organism ----------
  CELLTYPEPREDICT__res_pairwise <- pseudoBulk_testing(processed_data, 
                                                      organism_covariate=TRUE,
                                                      pairwise=TRUE,
													  pieces = 500,
      												  edgeR_obj = edgeR_obj,
														partition = partition)
  save(CELLTYPEPREDICT__res_pairwise, file = out)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE_PREDICT
  CELLTYPEPREDICT__res_organism_celltype <- pseudoBulk_testing(processed_data, 
                                                               organism_covariate=FALSE,
                                                               pairwise=TRUE, 
                                                               testing_against_internal_organism = TRUE,
                                                               testing_against = 'var_organism',
															   edgeR_obj = edgeR_obj,
														partition = partition)
  save(CELLTYPEPREDICT__res_organism_celltype, file = out)
  }
} else if (grepl('C', comp)){
  ##############################
  # now against cluster ---------
  ##############################
  if (grepl('w', comp)){
	piece = 1
  } else if (comp == 'C2') {
	piece = 500 
  } else {piece = 25} 
  if (grepl('1', comp)){
  # cluster against remaining, controlling for organism ------
  CLUSTER__res_againstAll <- pseudoBulk_testing(processed_data, 
                                                organism_covariate=TRUE,
                                                pairwise=FALSE,
                                                testing_against = 'cluster',
									            pieces =  piece,
											    edgeR_obj = edgeR_obj,
														partition = partition)
  save(CLUSTER__res_againstAll, file = out)
  } else if (grepl('2', comp)){
  # cluster against each cluster (pairwise), controlling for organism ----------
  CLUSTER__res_pairwise <- pseudoBulk_testing(processed_data, 
                                              organism_covariate=TRUE,
                                              pairwise=TRUE,
                                              testing_against = 'cluster',
											  pieces = piece,
											  edgeR_obj = edgeR_obj,
														partition = partition)
  save(CLUSTER__res_pairwise, file = out)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CLUSTER
  CLUSTER__res_organism_celltype <- pseudoBulk_testing(processed_data, 
                                                       organism_covariate=FALSE,
                                                       pairwise=TRUE, 
                                                       testing_against_internal_organism = TRUE,
                                                       testing_against = 'var_organism',
													   pieces = piece,
												       edgeR_obj = edgeR_obj,
														partition = partition)
  save(CLUSTER__res_organism_celltype, file = out)
  }
} else if (grepl('D', comp)){
  
  ##############################
  # SubCluster -------------
  # all testing only done within a cluster
  ##############################
  SUBCLUSTER__res_againstAll_list <- list()
  SUBCLUSTER__res_pairwise_list <- list()
  SUBCLUSTER__res_organism_celltype_list <- list()
  for (i in chunk(1:length(unique(umap$cluster)), 20)[[partition]]){
    umapT <- umap %>% filter(cluster == i)
	matT <- mat[,umapT$Barcode]
    info <- DataFrame(sample=as.factor(umapT$batch),
                      cluster = as.factor(umapT$subcluster),
                      organism = as.factor(umapT$organism))
    sum_mat5 <- sumCountsAcrossCells(matT, info, BPPARAM = multicoreParam)
    
    processed_data <- try({processing(sum_mat5, 
                                  testing_against = 'cluster') })

    # subcluster against remaining, controlling for organism ------
    if (grepl('1', comp) && class(processed_data) != 'try-error'){
    SUBCLUSTER__res_againstAll_list[[i]] <- try({pseudoBulk_testing(processed_data, 
                                                               organism_covariate=TRUE,
                                                               pairwise=FALSE,
                                                               testing_against = 'cluster',
											                   edgeR_obj = edgeR_obj,
														partition = FALSE) })
    } else if (grepl('2', comp) && class(processed_data) != 'try-error'){
    # subcluster against each subcluster (pairwise), controlling for organism ----------
    SUBCLUSTER__res_pairwise_list[[i]] <- try({pseudoBulk_testing(processed_data, 
                                                             organism_covariate=TRUE,
                                                             pairwise=TRUE,
                                                             testing_against = 'cluster',
														     edgeR_obj = edgeR_obj,
														partition = FALSE) })
    } else if (grepl('3', comp) && class(processed_data) != 'try-error'){
    # species against species, WITHIN A SUBCLUSTER
    SUBCLUSTER__res_organism_celltype_list[[i]] <- try({pseudoBulk_testing(processed_data, 
                                                                      organism_covariate=FALSE,
                                                                      pairwise=TRUE, 
                                                                      testing_against_internal_organism = TRUE,
                                                                      testing_against = 'var_organism',
																	  edgeR_obj = edgeR_obj,
														partition = FALSE) })
    }
  }
  if (grepl('1', comp)){
	for (i in seq_along(1:length(SUBCLUSTER__res_againstAll_list))){if (class(SUBCLUSTER__res_againstAll_list[[i]]) == 'try-error'){ SUBCLUSTER__res_againstAll_list[[i]] <- NULL  } }
    SUBCLUSTER__res_againstAll <- SUBCLUSTER__res_againstAll_list %>% bind_rows()
    save(SUBCLUSTER__res_againstAll, file = out)
  } else if (grepl('2', comp)){
	for (i in seq_along(1:length(SUBCLUSTER__res_pairwise_list))){if (class(SUBCLUSTER__res_pairwise_list[[i]]) == 'try-error'){ SUBCLUSTER__res_pairwise_list[[i]] <- NULL  } }
    SUBCLUSTER__res_pairwise <- SUBCLUSTER__res_pairwise_list %>% bind_rows()
	save(SUBCLUSTER__res_pairwise, file = out)
  } else if (grepl('3', comp)){ 
    for (i in seq_along(1:length(SUBCLUSTER__res_organism_celltype_list))){if (class(SUBCLUSTER__res_organism_celltype_list[[i]]) == 'try-error'){ SUBCLUSTER__res_organism_celltype_list[[i]] <- NULL  } } 
	SUBCLUSTER__res_organism_celltype <- SUBCLUSTER__res_organism_celltype_list %>% bind_rows()
    save(SUBCLUSTER__res_organism_celltype, out)
  }
}


