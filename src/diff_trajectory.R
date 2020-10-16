library(tidyverse)
library(scran)
library(TSCAN)
library(slingshot)
library(tictoc)
args = commandArgs(trailingOnly=TRUE)
load(args[1]) # slingshot obj

pt_tester <- function(sce, pseudoTime, block){
  pT <- colData(sce)[,pseudoTime] 
  block <- colData(sce)[,block]
  # remove NA pseudotimes
  pT_clean <- pT[!is.na(pT)]
  block_clean <- block[!is.na(pT)]
  # aggregate low n (<100) blocks into 1
  low_n_blocks <- block_clean %>% table() %>% enframe() %>% filter(value < 200) %>% pull(name)
  new_name <- paste(low_n_blocks, collapse = '__')
  block_clean[block_clean %in% low_n_blocks] <- new_name
  # diff test
  diff_PT <- TSCAN::testPseudotime(sce[,!is.na(pT)], 
                                   pseudotime = pT_clean, 
                                   block = block_clean,
								   BPPARAM = BiocParallel::MulticoreParam(12))
  diff_PT
}

sce <- sling$sling
diffPT <- list()
for (i in colnames(colData(sce)) %>% grep('slingP', .,  value = TRUE)){
  print(i)
  tic()
  diffPT[[i]] <- pt_tester(sce, i, 'study_accession')
  toc()
}

save(diffPT, file = args[2])
