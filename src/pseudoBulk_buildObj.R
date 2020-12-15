### 
# Use >= R/4.0
library(tictoc)
library(tidyverse)
library(scater)
library(Seurat)
library(edgeR)
library(scuttle)
library(BiocParallel)
multicoreParam <- MulticoreParam(workers = 12)

conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')

args <- commandArgs(trailingOnly = TRUE)

load(args[1]) #load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__onlyDROPLET__batch__scVI__dims6__preFilter__mindist0.1__nneighbors500.seuratObj.Rdata')
load(args[2])
load(args[3])
comp <- args[4]
#partition = args[4] %>% as.numeric()
out <- args[5]
###############
# functions -------
###############

source( paste0(git_dir, '/src/pseudoBulk_functions.R'))
#####################
# cluster into droplet umap
#####################
colnames(meta)[2] <- 'cluster'
umap <- umap %>% select(-cluster) %>% left_join(meta[,c(1,2)], by = 'Barcode')
####################
mat <- integrated_obj@assays$RNA@counts
mat <- mat[,umap$Barcode]
rm(integrated_obj)

umap <- umap %>% 
	#mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>%
	#mutate(CellType_predict = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType_predict)) %>%
	mutate(CTall = case_when(!is.na(CellType) ~ CellType, 
							!is.na(TabulaMurisCellType) ~ TabulaMurisCellType)) %>%
	mutate(CT_p_all = case_when(!is.na(CellType_predict) ~ CellType_predict, 
							!is.na(TabulaMurisCellType_predict) ~ TabulaMurisCellType)) %>%
	mutate(CTall = case_when(!grepl('Doub|Margin', CTall) ~ CTall))
if (grepl('A', comp)){
  ######################
  # celltype (pre-labelled/published) -------
  ####################
  info <- DataFrame(sample=as.factor(umap$batch),
                    celltype = as.factor(umap$CTall),
                    organism = as.factor(umap$organism))
  sum_mat2 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)
  
  #processed_data1 <- processing(sum_mat1)
  processed_data2 <- processing(sum_mat2)
  
  if (grepl('1', comp)){
  # celltype against remaining, controlling for organism ------
  edgeR_obj <- pseudoBulk_testing(processed_data2, 
                                                 organism_covariate=TRUE,
                                                 pairwise=FALSE,
												 save_edgeR_obj = TRUE)
  } else if (grepl('2', comp)){
  # celltype against each celltype (pairwise), controlling for organism ----------
  edgeR_obj <- pseudoBulk_testing(processed_data2, 
                                               organism_covariate=TRUE,
                                               pairwise=TRUE,
											   save_edgeR_obj = TRUE)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE
  edgeR_obj <- try({pseudoBulk_testing(processed_data2, 
                                                        organism_covariate=FALSE,
                                                        pairwise=TRUE, 
                                                        testing_against_internal_organism = TRUE,
                                                        testing_against = 'var_organism',
								                        save_edgeR_obj = TRUE) })
  }
} else if (grepl('B', comp)){
  ##############################
  # same, but with celltype predict (machine label all cells with label) ----------
  ##############################
  info <- DataFrame(sample=as.factor(umap$batch),
                    celltype = as.factor(umap$CT_p_all),
                    organism = as.factor(umap$organism))
  sum_mat3 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)
  
  processed_data3 <- processing(sum_mat3)
  if (grepl('1', comp)){
  # CellType_predict against remaining, controlling for organism ------
  edgeR_obj <- pseudoBulk_testing(processed_data3, 
                                                        organism_covariate=TRUE,
                                                        pairwise=FALSE,
												        save_edgeR_obj = TRUE)
  } else if (grepl('2', comp)){
  # CellType_predict against each celltype (pairwise), controlling for organism ----------
  edgeR_obj <- pseudoBulk_testing(processed_data3, 
                                                      organism_covariate=TRUE,
                                                      pairwise=TRUE,
      												  save_edgeR_obj = TRUE)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE_PREDICT
  edgeR_obj <- pseudoBulk_testing(processed_data3, 
                                                               organism_covariate=FALSE,
                                                               pairwise=TRUE, 
                                                               testing_against_internal_organism = TRUE,
                                                               testing_against = 'var_organism',
															   save_edgeR_obj = TRUE)
  }
} else if (grepl('C', comp)){
  ##############################
  # now against cluster ---------
  ##############################
  info <- DataFrame(sample=as.factor(umap$batch),
                    cluster = as.factor(umap$cluster),
                    organism = as.factor(umap$organism))
  sum_mat4 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)
  
  processed_data4 <- processing(sum_mat4, 
                                testing_against = 'cluster')

  piece = 500
  if (grepl('w', comp)){
	piece = 1
  }

  if (grepl('1', comp)){
  # cluster against remaining, controlling for organism ------
  edgeR_obj <- pseudoBulk_testing(processed_data4, 
                                                organism_covariate=TRUE,
                                                pairwise=FALSE,
                                                testing_against = 'cluster',
								                pieces = piece,
											    save_edgeR_obj = TRUE)
  } else if (grepl('2', comp)){
  # cluster against each cluster (pairwise), controlling for organism ----------
  edgeR_obj <- pseudoBulk_testing(processed_data4, 
                                              organism_covariate=TRUE,
                                              pairwise=TRUE,
                                              testing_against = 'cluster',
											  pieces = piece,
											  save_edgeR_obj = TRUE)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CLUSTER
  edgeR_obj <- pseudoBulk_testing(processed_data4, 
                                                       organism_covariate=FALSE,
                                                       pairwise=TRUE, 
                                                       testing_against_internal_organism = TRUE,
													   pieces = piece,
                                                       testing_against = 'var_organism',
												       save_edgeR_obj = TRUE)
  }
} else if (grepl('D', comp)){
  
  ##############################
  # SubCluster -------------
  # all testing only done within a cluster
  ##############################
  edgeR_obj <- list()
  for (i in chunk(1:length(unique(umap$cluster)), 20)[[partition]]){
    umapT <- umap %>% filter(cluster == i)
	matT <- mat[,umapT$Barcode]
    info <- DataFrame(sample=as.factor(umapT$batch),
                      cluster = as.factor(umapT$subcluster),
                      organism = as.factor(umapT$organism))
    sum_mat5 <- sumCountsAcrossCells(matT, info, BPPARAM = multicoreParam)
    
    processed_data5 <- try({processing(sum_mat5, 
                                  testing_against = 'cluster') })

    # subcluster against remaining, controlling for organism ------
    if (grepl('1', comp) && class(processed_data5) != 'try-error'){
    edgeR_obj[[i]] <- try({pseudoBulk_testing(processed_data5, 
                                                               organism_covariate=TRUE,
                                                               pairwise=FALSE,
                                                               testing_against = 'cluster',
											                   save_edgeR_obj = TRUE) })
    } else if (grepl('2', comp) && class(processed_data5) != 'try-error'){
    # subcluster against each subcluster (pairwise), controlling for organism ----------
    edgeR_obj[[i]] <- try({pseudoBulk_testing(processed_data5, 
                                                             organism_covariate=TRUE,
                                                             pairwise=TRUE,
                                                             testing_against = 'cluster',
														     save_edgeR_obj = TRUE) })
    } else if (grepl('3', comp) && class(processed_data5) != 'try-error'){
    # species against species, WITHIN A SUBCLUSTER
    edgeR_obj[[i]] <- try({pseudoBulk_testing(processed_data5, 
                                                                      organism_covariate=FALSE,
                                                                      pairwise=TRUE, 
                                                                      testing_against_internal_organism = TRUE,
                                                                      testing_against = 'var_organism',
																	  save_edgeR_obj = TRUE) })
    }
  }
}

save(edgeR_obj, file = out)
