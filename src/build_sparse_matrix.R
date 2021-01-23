# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)
#save(args, file = 'testing/build_sp_mtx_hm_nn.rdata')
library(Matrix)
library(tidyverse)
library(Seurat)
library(Matrix.utils)

species <- args[1]
cell_info_file <- args[2]
final_sparse_matrix_file <- args[3]
intron_sparse_matrix_file = args[4] 
metadata <- read_tsv(args[5])
gtf <-  rtracklayer::readGFF(args[6])
args_noidx = args[-(1:6)]
well_samples <-args_noidx[grepl('__counts.Rdata', args_noidx)]
spliced_droplet_samples <- args_noidx[grepl('matrix.Rdata', args_noidx)]
load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}


if (length(well_samples) > 0){

  all_well_data <- lapply(well_samples, load_rdata) %>% reduce(RowMergeSparseMatrices) 
  ## make row names for count (well) upper case
} else{
 all_well_data <- NULL 
}

# roll through UMI data, 
# correct gene names (upper case), force to be unique
# add sample ID to UMI
# e.g. AAATATAAAA_SRS2341234
empty_droplets <- list()
read_all_droplet_data = function(droplet_samples){
  sc_data <- list()
  droplet_sample_accessions <- list()

  i <- 1
  print(droplet_samples)
  for (sample in droplet_samples ){
	print(sample)
    drops <- load_rdata(sample)
    row.names(drops) <- row.names(drops) %>% str_remove_all('\\.$')
    colnames(drops) <- make.unique(colnames(drops))
    row.names(drops) <- make.unique(row.names(drops))
    sc_data[[i]] <- drops
    i <- i+1
    }
  
  all_droplet_data <- purrr::reduce(sc_data, RowMergeSparseMatrices)
  return(all_droplet_data)
}
# create naive fully merged  
all_spliced_droplet_data <- read_all_droplet_data(spliced_droplet_samples)
############################################################
# NOTE: Because sum(well_counts) >>> sum(droplet_counts), spliced/vs unspliced ratio will *look* off
# but its accurate. Double check 
################



## if well data exists, add it to droplet 
if (all(!is.null(all_well_data)) ){
  all_data <- RowMergeSparseMatrices(all_spliced_droplet_data, all_well_data)
} else {
  all_data <- all_spliced_droplet_data 
}

if(all(grepl('\\.\\d+$', sample(rownames(all_data), 10))) ){
  # if the gene ids have versions 
  rownames(all_data) <- str_remove_all(rownames(all_data), '\\.\\d+$')
}

all_data  <- all_data [row.names(all_data ) != 'fill.x', ] 
# create sample table
cell_info <- colnames(all_data ) %>% enframe() %>% 
  mutate(sample_accession = str_split(value, ':') %>% sapply(function(x) x[2]) %>% str_remove_all('\\.\\d+$'), 
         value = str_replace_all(value, ':','_')) %>% 
  left_join(metadata %>% select(-run_accession) %>% distinct()) %>% 
  data.frame()

row.names(cell_info) <- cell_info$value


cell_info <- cell_info %>% mutate(batch = paste(study_accession, Platform, Covariate, sep = '_'),
                                  batch2 = paste(study_accession, Covariate, sep = '_'),
                                  batch3 = paste(Platform, Covariate, sep = '_'))
#cell_info <- cell_info %>% mutate(Age = case_when(Age > 100 ~ 30, TRUE ~ Age))
# save barcodes for labelling with published cell type assignment 
write_tsv(cell_info, path = cell_info_file)

colnames(all_data) <- str_replace_all(colnames(all_data), ':', '_')
save(all_data, file = final_sparse_matrix_file, compress = FALSE)
## load intron quant 

intron_droplet_samples = spliced_droplet_samples %>% str_replace_all('matrix.Rdata', 'unspliced_matrix.Rdata')
all_intron_droplets = read_all_droplet_data(intron_droplet_samples)
if(all(grepl('\\.\\d+$', sample(rownames(all_intron_droplets), 10))) ){
  # if the gene ids have versions 
  rownames(all_intron_droplets) <- str_remove_all(rownames(all_intron_droplets), '\\.\\d+$')
}

colnames(all_intron_droplets) <- str_replace_all(colnames(all_intron_droplets), ':', '_')

save(all_intron_droplets, file = intron_sparse_matrix_file)




















