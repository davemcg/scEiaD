### 
# Use >= R/4.0
library(tidyverse)
library(scater)
library(Seurat)
library(edgeR)
library(BiocParallel)
multicoreParam <- MulticoreParam(workers = 8)
load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__onlyDROPLET__batch__scVI__dims6__preFilter__mindist0.1__nneighbors500.seuratObj.Rdata')
load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__onlyDROPLET__batch__scVI__dims6__preFilter__mindist0.1__nneighbors500.umap.Rdata')

###############
# functions -------
###############
processing <- function(sum_mat, testing_against = 'celltype'){
  colData <- sum_mat@colData
  cts <- sum_mat@assays@data$sum
  colnames(cts) <- colData %>% as_tibble() %>% select(-ncells) %>% unite(uname) %>% pull(uname)
  
  # remove sample/celltypes with low n
  discarded <- colData$ncells < 20
  colData <- colData[!discarded,]
  cts <- cts[,!discarded]
  
  
  # rename rod bipolar to bipolar
  colData[,'celltype'] <- gsub('Rod Bipolar Cells', 'Bipolar Cells', colData[,'celltype'] )
  
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
                               testing_against = 'celltype'){
  
  
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
    try({res <- glmQLFTest(fit, contrast = cont0)})
    
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
  results %>% bind_rows()
} 



mat <- integrated_obj@assays$RNA@counts

# # 1 mat by sample and celltype
# info <- DataFrame(sample=as.factor(umap$batch),
#                   celltype = as.factor(umap$CellType))
# sum_mat1 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)
# 2 mat by sample, celltype, and organism
info <- DataFrame(sample=as.factor(umap$batch),
                  celltype = as.factor(umap$CellType),
                  organism = as.factor(umap$organism))
sum_mat2 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)

#processed_data1 <- processing(sum_mat1)
processed_data2 <- processing(sum_mat2)

######################
# celltype (pre-labelled/published)
####################
# celltype against remaining, controlling for organism ------
res_againstAll <- pseudoBulk_testing(processed_data2, 
                                     organism_covariate=TRUE,
                                     pairwise=FALSE)
# celltype against each celltype (pairwise), controlling for organism ----------
res_pairwise <- pseudoBulk_testing(processed_data2, 
                                   organism_covariate=TRUE,
                                   pairwise=TRUE)
# species against species, WITHIN A CELLTYPE
res_organism_celltype <- pseudoBulk_testing(processed_data2, 
                                            organism_covariate=FALSE,
                                            pairwise=TRUE, 
                                            testing_against_internal_organism = TRUE,
                                            testing_against = 'var_organism')


##############################
# same, but with celltype predict (machine label all cells with label)
##############################
info <- DataFrame(sample=as.factor(umap$batch),
                  celltype = as.factor(umap$CellType_predict),
                  organism = as.factor(umap$organism))
sum_mat3 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)

processed_data3 <- processing(sum_mat3)
# CellType_predict against remaining, controlling for organism ------
res_againstAll <- pseudoBulk_testing(processed_data3, 
                                     organism_covariate=TRUE,
                                     pairwise=FALSE)
# CellType_predict against each celltype (pairwise), controlling for organism ----------
res_pairwise <- pseudoBulk_testing(processed_data3, 
                                   organism_covariate=TRUE,
                                   pairwise=TRUE)
# species against species, WITHIN A CELLTYPE_PREDICT
res_organism_celltype <- pseudoBulk_testing(processed_data3, 
                                            organism_covariate=FALSE,
                                            pairwise=TRUE, 
                                            testing_against_internal_organism = TRUE,
                                            testing_against = 'var_organism')


##############################
# now against cluster
##############################
info <- DataFrame(sample=as.factor(umap$batch),
                  cluster = as.factor(umap$cluster),
                  organism = as.factor(umap$organism))
sum_mat4 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)

processed_data4 <- processing(sum_mat4, 
                              testing_against = 'cluster')
# cluster against remaining, controlling for organism ------
res_againstAll <- pseudoBulk_testing(processed_data3, 
                                     organism_covariate=TRUE,
                                     pairwise=FALSE,
                                     testing_against = 'cluster')
# cluster against each celltype (pairwise), controlling for organism ----------
res_pairwise <- pseudoBulk_testing(processed_data3, 
                                   organism_covariate=TRUE,
                                   pairwise=TRUE,
                                   testing_against = 'cluster')
# species against species, WITHIN A CELLTYPE
res_organism_celltype <- pseudoBulk_testing(processed_data3, 
                                            organism_covariate=FALSE,
                                            pairwise=TRUE, 
                                            testing_against_internal_organism = TRUE,
                                            testing_against = 'var_organism')


##############################
# SubCluster
# all testing only done within a cluster
##############################
info <- DataFrame(sample=as.factor(umap$batch),
                  cluster = as.factor(umap$SubCluster),
                  organism = as.factor(umap$organism))
#for (cluster in umap)
sum_mat5 <- sumCountsAcrossCells(mat, info, BPPARAM = multicoreParam)

processed_data5 <- processing(sum_mat5, 
                              testing_against = 'cluster')
# cluster against remaining, controlling for organism ------
res_againstAll <- pseudoBulk_testing(processed_data3, 
                                     organism_covariate=TRUE,
                                     pairwise=FALSE,
                                     testing_against = 'cluster')
# cluster against each celltype (pairwise), controlling for organism ----------
res_pairwise <- pseudoBulk_testing(processed_data3, 
                                   organism_covariate=TRUE,
                                   pairwise=TRUE,
                                   testing_against = 'cluster')
# species against species, WITHIN A CELLTYPE
res_organism_celltype <- pseudoBulk_testing(processed_data3, 
                                            organism_covariate=FALSE,
                                            pairwise=TRUE, 
                                            testing_against_internal_organism = TRUE,
                                            testing_against = 'var_organism')



