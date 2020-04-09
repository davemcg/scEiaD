library(monocle3)
library(tidyverse)

#load('monocle_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.monocle.Rdata')
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims30__mindist0.1__nneighbors100.umap.Rdata')

args <- commandArgs(trailingOnly = TRUE)
model = args[3]

files <- Sys.glob(paste0('diff_testing/*',model,'.monocle_diff.Rdata')) 

gene_fits_all <- as_tibble()
count = 0
for (f in files){
  count = count + 1
  load(f)
  print(paste(count, f))
  gene_fits_all <- bind_rows(gene_fits_all, gene_fits)
}

coefficient_table <- coefficient_table(gene_fits_all %>% filter(status == 'OK')) %>% select(-model, -model_summary)
save(coefficient_table, file = args[1])
#save(gene_fits_all, file = args[2])
