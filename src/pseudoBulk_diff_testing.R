### 
# Use >= R/4.0
library(tidyverse)
library(scater)
library(Seurat)
library(edgeR)
library(BiocParallel)
multicoreParam <- MulticoreParam(workers = 12)

args <- commandArgs(trailingOnly = TRUE)

load(args[1]) #load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__onlyDROPLET__batch__scVI__dims6__preFilter__mindist0.1__nneighbors500.seuratObj.Rdata')
load(args[2]) # load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__onlyDROPLET__batch__scVI__dims6__preFilter__mindist0.1__nneighbors500.umap.Rdata')
comp <- args[3]
partition = args[4] %>% as.numeric()
out <- args[5]
###############
# functions -------
###############

chunk <- function(x,n) split(x, factor(sort(rank(x) %% n)))

processing <- function(sum_mat, testing_against = 'celltype'){
  colData <- sum_mat@colData
  cts <- sum_mat@assays@data$sum
  colnames(cts) <- colData %>% as_tibble() %>% select(-ncells) %>% unite(uname) %>% pull(uname)
  
  # remove sample/celltypes with low n
  discarded <- colData$ncells < 20
  colData <- colData[!discarded,]
  cts <- cts[,!discarded]
  
 
  if ('celltype' %in% colnames(colData)){ 
  	# rename rod bipolar to bipolar
  	colData[,'celltype'] <- gsub('Rod Bipolar Cells', 'Bipolar Cells', colData[,'celltype'] )
  }
  # remove celltypes/var of interest with 1 or fewer remaining replicates
  celltypes_to_retain <-  grep('Doubl|Margin', (colData[,testing_against] %>% table() )[(colData[,testing_against] %>% table() ) > 1] %>% names(), value = TRUE, invert = TRUE)
  cts <- cts[,(colData[,testing_against] %in% celltypes_to_retain)]
  colData <- colData[(colData[,testing_against] %in% celltypes_to_retain),]
  
  # add combo testing_against(e.g. celltype)_organism column
  colData$var_organism <- 
    paste0(colData[,testing_against], '_', colData[,'organism'])
  #keep <- edgeR::filterByExpr(cts, group = colData[,testing_against])
  #cts <- cts[keep,]
  
  colData[,testing_against] <- colData[,testing_against] %>% as.character() %>% as.factor()
  
  out <- list()
  out$cts <- cts
  out$colData <- colData
  out
}

# regex pattern to remove whitespace/ and / to make names safe for contrasts
# testing_against_internal_organism is for running the celltype_organism pairwise test
#   so you only test, for example amacrine_human against amacrine_mouse instead of
#   amacrine_human against rod_mouse
pseudoBulk_testing <- function(processed_data, 
                               organism_covariate = TRUE, 
                               pairwise = FALSE,
                               testing_against_internal_organism = FALSE,
                               regex_pattern = '\\s+|/', 
                               testing_against = 'celltype',
							   pieces = 20,
                               partition){
  
  
  # all pairwise combinations
  # remove white-space to facilitate contrast based results extraction
  if (pairwise){
    combinations <- combn(processed_data$colData[,testing_against] %>% unique() %>% gsub(regex_pattern,'', .), 2)
    combinations_fancy <- combn(processed_data$colData[,testing_against] %>% unique(), 2) # for naming
    # set up naming for testing for org specific genes changes
    # WITHIN a celltype (or cluster, etc)
    # e.g. cone human vs cone macaque
    if (testing_against_internal_organism) {
      comb_new <- data.frame(c(0,0))
      comb_fancy <- data.frame(c(0,0))
      for (i in seq_along(1:ncol(combinations))){
        # if celltype, etc in the first is the same as the second, then keep
        if (gsub('_.*','', combinations[,i][1]) == gsub('_.*','', combinations[,i][2])){
          comb_new <- cbind(comb_new, combinations[,i])
          comb_fancy <- cbind(comb_fancy, combinations_fancy[,i])
        }
      }
      comb_new <- comb_new[,-1]
      colnames(comb_new) <- seq_along(1:ncol(comb_new))
      combinations <- comb_new
      combinations_fancy <- comb_fancy[,-1]
    }
  } else {
    # if not pairwise, then test each var (e.g. cones in celltype) against the remaining (all other except cones)
    combinations <- processed_data$colData[,testing_against] %>% unique() %>% gsub(regex_pattern,'', .) %>% sort() %>% unique() %>% data.frame() %>% t()
  }
  
  # break into n pieces for parallelization faciliation on b2
  # take the n partition
  # skip for subcluster (partition set to FALSE)
  if (partition){
    combinations <- combinations[,chunk(1:ncol(combinations), pieces)[[partition]]]
    if (!pairwise){combinations <- combinations %>% t() %>% data.frame()}
  }
  # build factor for design matrix
  group <- as.factor(processed_data$colData[,testing_against] %>% gsub(regex_pattern,'', .))
  # study_covariate <- as.factor(colData_g1g2 %>% 
  #                                as_tibble() %>% 
  #                                separate(sample, into = c('study_accession', 'tech', 'batch'), remove = FALSE, sep = '_') %>% 
  #                                pull(study_accession))
  if (organism_covariate){
    org_covariate <- as.factor(processed_data$colData$organism %>% gsub(regex_pattern,'', .))
    design <- model.matrix(~0+group+org_covariate)
    colnames(design) <- gsub('^group','', colnames(design))
    colnames(design) <- gsub('^org_covariate','', colnames(design))
  } else {
    design <- model.matrix(~0+group)
    colnames(design) <- gsub('^group','', colnames(design))
  }
  
  y <- DGEList(counts = processed_data$cts, samples= processed_data$colData)
  y <- calcNormFactors(y)
  print('Estimate Dispersions')
  y <- estimateDisp(y, design)
  print('GLM Fitting')
  fit <- glmQLFit(y, design, robust=TRUE)
  #summary(fit$var.prior)
  #summary(fit$df.prior)
  results <- list()
  for (i in seq_along(1:ncol(combinations))){
    print(combinations[,i] %>% as.character())
    # contrast of all levels set to 0
    if (organism_covariate){
      cont0 <- integer((group %>% levels() %>% length()) +
                         (org_covariate %>% levels() %>% length()) -
                         1)
      names(cont0) <- c(group %>% levels(),
                        (org_covariate %>% levels())[-1])
    } else {
      cont0 <- integer(group %>% levels() %>% length())
      names(cont0) <- group %>% levels()
    }
    # now set contrasts for the pairwise test
    if (pairwise){
      cont0[combinations[,i][1]] = 1
      cont0[combinations[,i][2]] = -1
    } else { # and the one vs all test
      cont0[1:length(cont0)] = -(1/(length(levels(group))-1))
      cont0[combinations[,i]] = 1
      if (organism_covariate){
        cont0['Macacafascicularis'] = 0
        cont0['Musmusculus'] = 0
      }
    }
    res <- try({glmQLFTest(fit, contrast = cont0)})
    if (class(res) == 'try-error') {
      if (pairwise){test = paste0(combinations_fancy[,i][1], ' vs ', combinations_fancy[,i][2])} else {test = combinations[,i]}
      comparison = 'failure'
      results <- data.frame(cbind(test, failure)) %>% as_tibble()
    } else { 
      if (pairwise){
        results[[i]] <- topTags(res, n = 10000000) %>% 
          data.frame() %>% 
          as_tibble(rownames = 'Gene') %>% 
          mutate(test = paste0(combinations_fancy[,i][1], ' vs ', combinations_fancy[,i][2])) %>% 
          mutate(comparison = res$comparison)
      } else {
        results[[i]] <- topTags(res, n = 10000000) %>% 
          data.frame() %>% 
          as_tibble(rownames = 'Gene') %>% 
          mutate(test = combinations[,i]) %>% 
          mutate(comparison = res$comparison)
      }
   }
  }
  results %>% bind_rows()
} 



mat <- integrated_obj@assays$RNA@counts
mat <- mat[,umap$Barcode]
rm(integrated_obj)

if (grepl('A', comp)){
  ######################
  # celltype (pre-labelled/published) -------
  ####################
  info <- DataFrame(sample=as.factor(umap$batch),
                    celltype = as.factor(umap$CellType),
                    organism = as.factor(umap$organism))
  sum_mat2 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)
  
  #processed_data1 <- processing(sum_mat1)
  processed_data2 <- processing(sum_mat2)
  
  if (grepl('1', comp)){
  # celltype against remaining, controlling for organism ------
  CELLTYPE__res_againstAll <- pseudoBulk_testing(processed_data2, 
                                                 organism_covariate=TRUE,
                                                 pairwise=FALSE,
												 partition = partition)
  save(CELLTYPE__res_againstAll, file = out)
  } else if (grepl('2', comp)){
  # celltype against each celltype (pairwise), controlling for organism ----------
  CELLTYPE__res_pairwise <- pseudoBulk_testing(processed_data2, 
                                               organism_covariate=TRUE,
                                               pairwise=TRUE,
											   partition = partition)
  save(CELLTYPE__res_pairwise, file = out)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE
  CELLTYPE__res_organism_celltype <- pseudoBulk_testing(processed_data2, 
                                                        organism_covariate=FALSE,
                                                        pairwise=TRUE, 
                                                        testing_against_internal_organism = TRUE,
                                                        testing_against = 'var_organism',
								                        partition = partition)
  save(CELLTYPE__res_organism_celltype, file = out)
  }
} else if (grepl('B', comp)){
  ##############################
  # same, but with celltype predict (machine label all cells with label) ----------
  ##############################
  info <- DataFrame(sample=as.factor(umap$batch),
                    celltype = as.factor(umap$CellType_predict),
                    organism = as.factor(umap$organism))
  sum_mat3 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)
  
  processed_data3 <- processing(sum_mat3)
  if (grepl('1', comp)){
  # CellType_predict against remaining, controlling for organism ------
  CELLTYPEPREDICT__res_againstAll <- pseudoBulk_testing(processed_data3, 
                                                        organism_covariate=TRUE,
                                                        pairwise=FALSE,
												        partition = partition)
  save(CELLTYPEPREDICT__res_againstAll, file = out)
  } else if (grepl('2', comp)){
  # CellType_predict against each celltype (pairwise), controlling for organism ----------
  CELLTYPEPREDICT__res_pairwise <- pseudoBulk_testing(processed_data3, 
                                                      organism_covariate=TRUE,
                                                      pairwise=TRUE,
      												  partition = partition)
  save(CELLTYPEPREDICT__res_pairwise, file = out)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE_PREDICT
  CELLTYPEPREDICT__res_organism_celltype <- pseudoBulk_testing(processed_data3, 
                                                               organism_covariate=FALSE,
                                                               pairwise=TRUE, 
                                                               testing_against_internal_organism = TRUE,
                                                               testing_against = 'var_organism',
															   partition = partition)
  save(CELLTYPEPREDICT__res_organism_celltype, file = out)
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
  if (grepl('1', comp)){
  # cluster against remaining, controlling for organism ------
  CLUSTER__res_againstAll <- pseudoBulk_testing(processed_data4, 
                                                organism_covariate=TRUE,
                                                pairwise=FALSE,
                                                testing_against = 'cluster',
											    partition = partition)
  save(CLUSTER__res_againstAll, file = out)
  } else if (grepl('2', comp)){
  # cluster against each celltype (pairwise), controlling for organism ----------
  CLUSTER__res_pairwise <- pseudoBulk_testing(processed_data4, 
                                              organism_covariate=TRUE,
                                              pairwise=TRUE,
                                              testing_against = 'cluster',
											  partition = partition)
  save(CLUSTER__res_pairwise, file = out)
  } else if (grepl('3', comp)){
  # species against species, WITHIN A CELLTYPE
  CLUSTER__res_organism_celltype <- pseudoBulk_testing(processed_data4, 
                                                       organism_covariate=FALSE,
                                                       pairwise=TRUE, 
                                                       testing_against_internal_organism = TRUE,
                                                       testing_against = 'var_organism',
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
    
    processed_data5 <- try({processing(sum_mat5, 
                                  testing_against = 'cluster') })

    # subcluster against remaining, controlling for organism ------
    if (grepl('1', comp) && class(processed_data5) != 'try-error'){
    SUBCLUSTER__res_againstAll_list[[i]] <- try({pseudoBulk_testing(processed_data5, 
                                                               organism_covariate=TRUE,
                                                               pairwise=FALSE,
                                                               testing_against = 'cluster',
											                   partition = FALSE) })
    } else if (grepl('2', comp) && class(processed_data5) != 'try-error'){
    # subcluster against each subcluster (pairwise), controlling for organism ----------
    SUBCLUSTER__res_pairwise_list[[i]] <- try({pseudoBulk_testing(processed_data5, 
                                                             organism_covariate=TRUE,
                                                             pairwise=TRUE,
                                                             testing_against = 'cluster',
														     partition = FALSE) })
    } else if (grepl('3', comp) && class(processed_data5) != 'try-error'){
    # species against species, WITHIN A SUBCLUSTER
    SUBCLUSTER__res_organism_celltype_list[[i]] <- try({pseudoBulk_testing(processed_data5, 
                                                                      organism_covariate=FALSE,
                                                                      pairwise=TRUE, 
                                                                      testing_against_internal_organism = TRUE,
                                                                      testing_against = 'var_organism',
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


