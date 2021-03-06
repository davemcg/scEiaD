library(monocle3)
library(tidyverse)

#load('monocle_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.monocle.Rdata')
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.umap.Rdata')

args <- commandArgs(trailingOnly = TRUE)
#load(args[1]) # seurat cluster obj
load(args[1]) # monocle obj
#cds_retina@colData$monocleCluster <- cds_retina@clusters$UMAP$clusters %>% as.factor()
#cds_retina@colData$seuratCluster <- meta[,2] %>% pull(1) %>%  as.factor()
#cds_retina@colData$seuratSubCluster <- meta[,3] %>% pull(1) %>%  as.factor()

piece_n <- as.numeric(args[2]) # number of pieces to subset into
the_n <- as.numeric(args[3]) # which subset this run is doing
model <- args[4] #e.g "cluster-batch-percent.mt" where cluster is the comparison of interest and the remainder are covariates
if (model == 'A'){
  model <- '~seuratCluster+batch+percent.mt'
} else if (model == 'B') {
  model <- '~seuratCluster+batch+percent.mt+organism'
} else if (model == 'C') {
  model <- '~CellType_predict+batch+percent.mt'
} else if (model == 'D') {
  model <- '~CellType_predict+batch+percent.mt+organism'
} else if (model == 'E') {
  model <- '~CellType+batch+percent.mt'
} else if (model == 'F') {
  model <- '~CellType+batch+percent.mt+organism'
} else if (model == 'G') {
  model <- '~seuratSubCluster+batch+percent.mt'
}
test_var <- gsub('~','',model) %>% gsub('\\+.*','', .)
output <- args[5]

chunks <- split(1:nrow(cds_retina), 
                sort(rep_len(1:piece_n, 
                             length(1:nrow(cds_retina)))
                     )
                ) 
print(chunks[[the_n]])

cds_subset <- cds_retina[chunks[[the_n]], !is.na(colData(cds_retina)[,test_var])]
if (test_var == 'seuratSubCluster'){
  gene_fits_list <- list()
  for (i in colData(cds_retina)$seuratCluster %>% as.character() %>% unique()){
    print(i)
    subsub <- cds_subset[, colData(cds_retina)$seuratCluster == i]
    fits <- fit_models(subsub, model_formula_str = model)
    gene_fits_list[[i]] <- fits
  }
  gene_fits <- bind_rows(gene_fits_list)
} else { gene_fits <- fit_models(cds_subset, model_formula_str = model)
}

save(gene_fits, file = output)
