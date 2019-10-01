# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(tidyverse)
library(Seurat)
species <- str_split(args[1], '/')[[1]][2]

metadata <- read_tsv(args[4])
tx <- read_tsv(args[5], col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

# well data
if (species != "Macaca_fascicularis"){
  load(args[5]) # well data
  rdata_files = args[6:length(args)]
  ## make row names for count (well) upper case
  row.names(count) <- toupper(row.names(count))
} else {
  rdata_files = args[5:length(args)]
}

# roll through UMI data, 
# correct gene names (upper case), force to be unique
# add sample ID to UMI
# e.g. AAATATAAAA_SRS2341234
sc_data <- list()
droplet_samples <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  droplet_samples <- c(droplet_samples, sample_accession)
  load(file)
  row.names(res_matrix) <- row.names(res_matrix) %>% 
    enframe(value = 'id') %>% 
    left_join(., tx, by = 'id') %>% 
    pull(gene) %>% 
    toupper()
  # remove cells which have more than 6000 quantified genes (likely doublets)
  #res_matrix <- res_matrix[,diff(res_matrix@p) < 6000]
  colnames(res_matrix) <- make.unique(colnames(res_matrix))
  colnames(res_matrix) <- paste0(colnames(res_matrix), "_", sample_accession)
  row.names(res_matrix) <- make.unique(row.names(res_matrix))
  sc_data[[sample_accession]] <- res_matrix
}

# create naive fully merged  
droplet <- Reduce(cbind, sc_data) 
if (species != "Macaca_fascicularis"){
  m <- Matrix.utils::merge.Matrix(count, droplet, by.x=row.names(count), by.y = row.names(droplet))
} else {
  m <- droplet 
}

m <- m[row.names(m) != 'fill.x', ] 
# create sample table
cell_info <- colnames(m) %>% enframe() %>% 
  mutate(sample_accession = str_extract(value, 'SRS\\d+')) %>% 
  left_join(metadata %>% select(-run_accession) %>% unique()) %>% 
  data.frame()
row.names(cell_info) <- cell_info$value

cell_info <- cell_info %>% mutate(batch = paste(study_accession, Platform, Covariate, sep = '_'),
                                  batch2 = paste(study_accession, Covariate, sep = '_'),
                                  batch3 = paste(Platform, Covariate, sep = '_'))
#cell_info <- cell_info %>% mutate(Age = case_when(Age > 100 ~ 30, TRUE ~ Age))
# save barcodes for labelling with published cell type assignment 
save(cell_info, file = args[1])
save(m, file = args[2], compress = FALSE)