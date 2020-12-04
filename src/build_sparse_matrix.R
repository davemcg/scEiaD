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
metadata <- read_tsv(args[4])
gtf <-  rtracklayer::readGFF(args[5])
args_noidx = args[-(1:5)]
well_samples <-args_noidx[grepl('__counts.Rdata', args_noidx)]
droplet_samples <- args_noidx[grepl('matrix.Rdata', args_noidx)]
#tx <- read_tsv(args[5], col_names = FALSE) %>% select(2,3) %>% unique()
#colnames(tx) <- c('id', 'gene')


load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  print(x)
  print(var)
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
sc_data <- list()
droplet_sample_accessions <- list()
empty_droplets <- list()

for (sample in droplet_samples ){
  drops <- load_rdata(sample)
  sample_accession = str_extract(sample, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')
  if (is.null(dim(drops))){
    empty_droplets <- c(empty_droplets, sample_accession)
  } else if (ncol(drops) == 0) {
    empty_droplets <- c(empty_droplets, sample_accession)
  } else {
  
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

# create naive fully merged  
all_droplet_data <- reduce(sc_data, RowMergeSparseMatrices)
if (species != "Macaca_fascicularis"){
  all_data <- RowMergeSparseMatrices(all_droplet_data, all_well_data)
} else {
  all_data <- all_droplet_data 
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
####VS:
# when generating alignment indices with BUSparse, geneid's are used; for well data, I could summarrise to gene name when all nonUMI is merged
# but felt that it was better to treat all data the same. multiple geneids map to the same gene_name so aggregate by gene name 
####
# geneid2gene_name <- gtf %>% filter(type == 'transcript') %>% select(gene_id, gene_name) %>% distinct
# gene_name_ordered_group <- row.names(all_data) %>% 
#   {tibble(gene_id = .)} %>% 
#   inner_join(geneid2gene_name) %>% 
#   pull(gene_name) %>% toupper()
# all_data_by_genename <- aggregate.Matrix(x = all_data, groupings = gene_name_ordered_group, fun = 'sum')

#cell_info <- cell_info %>% mutate(Age = case_when(Age > 100 ~ 30, TRUE ~ Age))
# save barcodes for labelling with published cell type assignment 
write_tsv(cell_info, path = cell_info_file)
save(all_data, file = final_sparse_matrix_file, compress = FALSE)
















