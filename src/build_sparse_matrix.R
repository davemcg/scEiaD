# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)

#args <- c('seurat_obj/Mus_musculus__standard_and_SCT__late__batch.seuratV3.Rdata','late','batch','/home/mcgaugheyd/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv','references/gencode.vM22.metadata.MGI_tx_mapping.tsv','quant/Mus_musculus/counts.Rdata','quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

library(Matrix)
library(tidyverse)
library(Seurat)

metadata <- read_tsv(args[3])
tx <- read_tsv(args[4], col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

# well data
load(args[5])
# # remove cells with > 10000 or < 1000 detected genes
# count <- count[,(diff(count@p) < 10000)]
# count <- count[,(diff(count@p) > 1000)]

# extract species
species <- str_split(args[6], '/')[[1]][2]

# sparse matrix files
rdata_files = args[8:length(args)]

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

# create naive fully merged for 
## make row names for count (well) upper case
row.names(count) <- toupper(row.names(count))
droplet <- Reduce(cbind, sc_data) 
m <- Matrix.utils::merge.Matrix(count, droplet, by.x=row.names(count), by.y = row.names(droplet))
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