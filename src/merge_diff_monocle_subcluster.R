library(tidyverse)

#load('monocle_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.monocle.Rdata')
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.umap.Rdata')

args <- commandArgs(trailingOnly = TRUE)
model = args[2]

#files <- Sys.glob(paste0('diff_testing/*',model,'.diff.Rdata'))
files <- args[3:length(args)]

all = list()
for (f in files){
  load(f)
  all[[f]] <- coefficient_table
}

coefficient_table <- all %>% bind_rows()
save(coefficient_table, file = args[1])
#save(gene_fits_all, file = args[2])
~
