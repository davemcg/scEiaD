### 
# Use >= R/4.0
library(tictoc)
library(tidyverse)
library(scater)
library(Seurat)
library(edgeR)
library(BiocParallel)
multicoreParam <- MulticoreParam(workers = 12)

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
							   pieces = 25,
                               partition = NA,
							   edgeR_obj = NA,
							   save_edgeR_obj = FALSE){
  
  ocular_celltypes <- c('AC/HC_Precurs','Amacrine Cells','Artery','Astrocytes','B-Cell','Bipolar Cells','Choriocapillaris','Cones','Early RPCs','Endothelial','Fibroblasts','Horizontal Cells','Late RPCs','Macrophage','Mast','Melanocytes','Microglia','Muller Glia','Neurogenic Cells','Pericytes','Photoreceptor Precursors','Red Blood Cells','Retinal Ganglion Cells','Rod Bipolar Cells','Rods','RPCs','RPE','Schwann','Smooth Muscle Cell','T-Cell','Vein') 
  # all pairwise combinations
  # remove white-space to facilitate contrast based results extraction
  if (pairwise){
    combinations <- combn(processed_data$colData[,testing_against] %>% unique() %>% gsub(regex_pattern,'', .), 2)
    combinations_fancy <- combn(processed_data$colData[,testing_against] %>% unique() %>% as.character(), 2) # for naming
    
	# remove non-ocular vs non-ocular tests
	if (testing_against == 'celltype'){
		c2 <- combinations[,apply(combinations, 2, function(x) sum(x %in% ocular_celltypes)) > 0 ]
		c2_fancy <- combinations_fancy[,apply(combinations, 2, function(x) sum(x %in% ocular_celltypes)) > 0 ]
		combinations <- c2; combinations_fancy <- c2_fancy
	}
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
  	combinations_fancy <- combinations
  }
  
  # break into n pieces for parallelization faciliation on b2
  # take the n partition
  # skip for subcluster (partition set to FALSE)
  if (!is.na(partition)){
    combinations <- combinations[,chunk(1:ncol(combinations), pieces)[[partition]]] %>% data.frame()
	combinations_fancy <- combinations_fancy[,chunk(1:ncol(combinations_fancy), pieces)[[partition]]] %>% data.frame()
    if (!pairwise){
		combinations <- combinations %>% t() %>% data.frame()
		combinations_fancy <- combinations_fancy %>% t() %>% data.frame()
	}
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
  
  if (is.na(edgeR_obj)) {
    y <- DGEList(counts = processed_data$cts, samples= processed_data$colData)
    y <- calcNormFactors(y)
    print('Estimate Dispersions')
    y <- estimateDisp(y, design)
    print('GLM Fitting')
    fit <- glmQLFit(y, design, robust=TRUE)
  } else {
	print("Using provided edgeR obj!")
	y <- edgeR_obj$y
	fit <- edgeR_obj$fit
} 

  if (save_edgeR_obj) {
	print("Saving edgeR obj!")
	edgeR <- list()
	edgeR$y <- y
	edgeR$fit <- fit
	edgeR$processed_data <- processed_data 
	edgeR
  } else {
	print("Diff testing begins!")
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
        if (organism_covariate &
				sum(grepl('Maca', processed_data$colData$var_organism %>% 
					unique()) %>% sum() ) != 0){
          cont0['Macacafascicularis'] = 0
          cont0['Musmusculus'] = 0
        } else {
		  cont0['Musmusculus'] = 0
		}
      }
      tic(); res <- try({glmQLFTest(fit, contrast = cont0)}); toc()
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
} 

