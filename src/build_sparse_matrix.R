# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)
#save(args, file = 'testing/build_sp_mtx_hm.rdata')
library(Matrix)
library(tidyverse)
library(Seurat)
library(Matrix)
library(Matrix.utils)

species <- args[1]
cell_info_file <- args[2]
final_sparse_matrix_file <- args[3]
intron_sparse_matrix_file = args[4] 
empty_droplet_rdata = args[5]
metadata <- read_tsv(args[6])
gtf <-  rtracklayer::readGFF(args[7])
args_noidx = args[-(1:7)]
well_samples <-args_noidx[grepl('__counts.Rdata', args_noidx)]
spliced_droplet_samples <- args_noidx[grepl('matrix.Rdata', args_noidx)]
#tx <- read_tsv(args[5], col_names = FALSE) %>% select(2,3) %>% unique()
#colnames(tx) <- c('id', 'gene')


load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}


print(species)
if (species != "Macaca_fascicularis"){

  all_well_data <- lapply(well_samples, load_rdata) %>% reduce(RowMergeSparseMatrices) 
  ## make row names for count (well) upper case
} 

# roll through UMI data, 
# correct gene names (upper case), force to be unique
# add sample ID to UMI
# e.g. AAATATAAAA_SRS2341234
empty_droplets <- list()
read_all_droplet_data = function(droplet_samples){
  sc_data <- list()
  droplet_sample_accessions <- list()


  for (sample in droplet_samples ){
    drops <- load_rdata(sample)
    sample_accession = str_extract(sample, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')
    if (is.null(dim(drops))){
      empty_droplets <- c(empty_droplets, sample_accession)# this is coming in from global env
    } else{
    
    droplet_sample_accessions<- c(droplet_sample_accessions, sample_accession)

    #if (species != "Macaca_fascicularis") {
    row.names(drops) <- row.names(drops) %>% str_remove_all('\\.$')
      
    # } 
    colnames(drops) <- make.unique(colnames(drops))
    colnames(drops) <- paste0(colnames(drops), "_", sample_accession)
    row.names(drops) <- make.unique(row.names(drops))
    sc_data[[sample_accession]] <- drops
    }
  }
  all_droplet_data <- reduce(sc_data, RowMergeSparseMatrices)
  return(all_droplet_data)
}
# create naive fully merged  
all_spliced_droplet_data <- read_all_droplet_data(spliced_droplet_samples)
if (species != "Macaca_fascicularis"){
  all_data <- RowMergeSparseMatrices(all_spliced_droplet_data, all_well_data)
} else {
  all_data <- all_spliced_droplet_data 
  if (!grepl('hs-homo_sapiens', args[5])){
    row.names(all_data) <- row.names(all_data) %>% str_remove('\\.\\d+$')
    # macaque gtf doesnt have gene versions 
  }
}

all_data  <- all_data [row.names(all_data ) != 'fill.x', ] 
# create sample table
cell_info <- colnames(all_data ) %>% enframe() %>% 
  mutate(sample_accession = str_extract(value, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')) %>% 
  left_join(metadata %>% select(-run_accession) %>% unique()) %>% 
  data.frame()
row.names(cell_info) <- cell_info$value

cell_info <- cell_info %>% mutate(batch = paste(study_accession, Platform, Covariate, sep = '_'),
                                  batch2 = paste(study_accession, Covariate, sep = '_'),
                                  batch3 = paste(Platform, Covariate, sep = '_'))
#cell_info <- cell_info %>% mutate(Age = case_when(Age > 100 ~ 30, TRUE ~ Age))
# save barcodes for labelling with published cell type assignment 
write_tsv(cell_info, path = cell_info_file)
save(all_data, file = final_sparse_matrix_file, compress = FALSE)
## load intron quant 
all_empty_droplets = list(spliced = empty_droplets)
empty_droplets = list()

intron_droplet_samples = spliced_droplet_samples %>% str_replace_all('matrix.Rdata', 'unspliced_matrix.Rdata')
all_intron_droplets = read_all_droplet_data(intron_droplet_samples)
if (species == "Macaca_fascicularis"){
   
  if (!grepl('hs-homo_sapiens', args[5])){
    row.names(all_intron_droplets) <- row.names(all_intron_droplets) %>% str_remove('\\.\\d+$')
    # macaque gtf doesnt have gene versions 
  }
}
save(all_intron_droplets, file = intron_sparse_matrix_file)
all_empty_droplets[['intron']] = empty_droplets
save(all_empty_droplets, file = empty_droplet_rdata )




















