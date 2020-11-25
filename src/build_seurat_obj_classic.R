# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)
save(args, file ='testing/bso_args.Rdata')
#Sys.setenv(SCIAD_CONFIG = '/data/swamyvs/scEiaD/config.yaml')
#args <- c('seurat_obj/Mus_musculus__standard_and_SCT__late__batch.seuratV3.Rdata','late','batch','/home/mcgaugheyd/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv','references/gencode.vM22.metadata.MGI_tx_mapping.tsv','quant/Mus_musculus/counts.Rdata','quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

library(Matrix)
library(Matrix.utils)
library(tidyverse)
library(Seurat)
library(scran)
library(future)
library(quminorm)
library(glue)
library(yaml)
plan(strategy = "multicore", workers = 4)
config=read_yaml(Sys.getenv('SCIAD_CONFIG'))
git_dir=config$git_dir
working_dir=config$working_dir
setwd(working_dir)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)
print('load args')
set = args[2] # early, late, full, downsampled
covariate = args[3] # study_accession, batch, etc.
transform = args[4] # SCT or standard seurat
combination = args[5] # mouse, mouse and macaque, mouse and macaque and human
n_features = args[6] %>% as.numeric()
cell_info <- read_tsv(args[7]) # cell_info.tsv
cell_info$batch <- gsub(' ', '', cell_info$batch)
# set batch covariate for well data to NA, as any splits risks making the set too small
print('cell info import')
cell_info <- cell_info %>% 
  mutate(batch = case_when(study_accession == 'SRP125998' ~ paste0(study_accession, "_", Platform, '_NA'),
                           TRUE ~ batch)) %>% 
  mutate(batch = gsub(' ', '_', batch))

load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}

m <- load_rdata(args[8])

# n=100
# 
# 
# if(n >0){
#     top_n_genes_hs <- rdata[['Homo_sapiens']] %>% 
#       {tibble(gene = row.names(.), count = rowSums(.) )} %>% 
#       filter(!gene%in% shared_genes) %>% 
#       arrange(desc(count)) %>% 
#       head(n) %>% pull(gene)
#     top_n_genes_mm <- rdata[['Mus_musculus']] %>% 
#       {tibble(gene = row.names(.), count = rowSums(.) )} %>% 
#       filter(!gene%in% shared_genes) %>% 
#       arrange(desc(count)) %>% 
#       head(n) %>% pull(gene)
#     top_n_genes_mf <- rdata[['Macaca_fascicularis']] %>% 
#       {tibble(gene = row.names(.), count = rowSums(.) )} %>% 
#       filter(!gene%in% shared_genes) %>% 
#       arrange(desc(count)) %>% 
#       head(n) %>% pull(gene)
#     var_genes <- list('Homo_sapiens' = top_n_genes_hs, 
#                       'Mus_musculus' = top_n_genes_mm,
#                       'Macaca_fascicularis' = top_n_genes_mf
#                       )
# 
#     for (i in names(rdata)){
#       file_cut_down[[i]] <- rdata[[i]][c(shared_genes, var_genes[[i]]) ,]
#     }
#     m <- file_cut_down %>% purrr::reduce(RowMergeSparseMatrices)
#     
# } else{
# 
#   for (i in names(rdata)){
#     file_cut_down[[i]] <- rdata[[i]][shared_genes,]
#     
#   }
#   m <- file_cut_down %>% purrr::reduce(cbind)
#   
# }

# file_cut_down <- list()
# mito_list <- list()
# for (i in names(rdata)){
#   file_cut_down[[i]] <- rdata[[i]][shared_genes,]
# }
# m <- file_cut_down %>% purrr::reduce(cbind)

print('Splitting time')
# custom combos / sets
m_early <- m[,cell_info %>% filter(Age < 10) %>% pull(value)]
m_late <- m[,cell_info %>% filter(Age >= 10) %>% pull(value)]
m_test <- m[,sample(1:ncol(m), 10000)]


precursors <- c('AC/HC_Precurs', 'Early RPCs', 'Late RPCs', 'Neurogenic Cells', 'Photoreceptor Precursors', 'RPCs')
# load(config$mso_umap_file)
# if (set == 'cones') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Cones', precursors)) %>% pull(Barcode)]
# }
# if (set == 'bipolar') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Bipolar Cells', precursors)) %>% pull(Barcode)]
# }
# if (set == 'rods') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Rods', precursors)) %>% pull(Barcode)]
# }
# if (set == 'mullerglia') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Muller Glia', precursors)) %>% pull(Barcode)]
# }
# if (set == 'amacrine') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Amacrine Cells', precursors)) %>% pull(Barcode)]
# }
# if (set == 'rgc') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Retinal Ganglion Cells', precursors)) %>% pull(Barcode)]
# }
# if (set == 'hc') {
# 	m_subset = m[, umap %>% filter(CellType_predict %in% c('Horizontal Cells', precursors)) %>% pull(Barcode)]
# }


m_onlyDROPLET <-  m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]
m_TABULA_DROPLET <- m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
m_onlyWELL <- m[,cell_info %>% filter(!Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]

m_onlyDROPLET_hs <-  m[,cell_info %>% filter(organism == 'Homo sapiens', Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]
m_onlyDROPLET_mm  <-  m[,cell_info %>% filter(organism == 'Mus musculus', Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]

downsample_samples <- 
  cell_info %>% 
  group_by(batch) %>% 
  sample_n(2000, replace = TRUE) %>% 
  unique() %>% 
  pull(value)
m_downsample <- m[,downsample_samples]

source(glue('{git_dir}/src/make_seurat_obj_functions.R') )

if (set == 'early'){
  print("Running Early")
  seurat__standard <- make_seurat_obj(m_early, split.by = covariate)
} else if (set == 'late'){
  print("Running Late")
  seurat__standard <- make_seurat_obj(m_late, split.by = covariate)
} else if (set == 'full'){
  print("Running Full")
  seurat__standard <- make_seurat_obj(m, split.by = covariate)
} else if (set == 'onlyDROPLET'){
  print("Running onlyDROPLET (remove well based)")
  seurat__standard <- make_seurat_obj(m_onlyDROPLET, split.by = covariate, keep_well = FALSE)
}  else if (set == 'TabulaDroplet'){
  print("Running onlyDROPLET with Tabula Muris (no well)")
  seurat__standard <- make_seurat_obj(m_TABULA_DROPLET, split.by = covariate, keep_well = FALSE)
} else if (set == 'HSdroplet') {
  seurat__standard <- make_seurat_obj(m_onlyDROPLET_hs, split.by = covariate, keep_well = FALSE)  
} else if (set == 'MMdroplet') {
  seurat__standard <- make_seurat_obj(m_onlyDROPLET_mm, split.by = covariate, keep_well = FALSE)
} else if (set == 'onlyWELL' & transform == 'counts'){
  print("Running onlyWELL with gene length correction as per scVI (remove droplet based)") 
  seurat__standard <- make_seurat_obj(m_onlyWELL, split.by = covariate, keep_droplet = FALSE, lengthCor = TRUE)
} else if (set == 'onlyWELL') {
  print("Running onlyWELL (remove droplet based)") 
  seurat__standard <- make_seurat_obj(m_onlyWELL, split.by = covariate, keep_droplet = FALSE)
} else if (set == 'downsample'){
  print("Running downsample")
  seurat__standard <- make_seurat_obj(m_downsample, split.by = covariate)
} else if (set == 'universe'){
  seurat__standard <- make_seurat_obj(m, split.by = covariate, lengthCor = TRUE)
} else if (set %in% c('cones', 'hc', 'rgc', 'amacrine', 'mullerglia', 'bipolar', 'rods' )){
  seurat__standard <- make_seurat_obj(m_subset, split.by = covariate, keep_well = FALSE)
} else if (set == 'raw') {
  seurat__standard <- m
}

if (transform == 'SCT'){
  s_data_list<- SplitObject(seurat__standard, split.by = covariate)
  seurat__SCT <- seurat_sct(s_data_list)
  save(seurat__SCT, 
       file = args[1], compress = FALSE)
} else if (transform == 'scran'){
  seurat__standard <- scran_norm(seurat__standard, split.by = covariate )
  save(seurat__standard, 
       file = args[1], compress = FALSE)
} else if (transform == 'libSize'){
  seurat__standard <- library.size.normalize(seurat__standard)
  save(seurat__standard, 
       file = args[1], compress = FALSE)
} else if (transform == 'sqrt') {
  seurat__standard <- library.size.normalize(seurat__standard, sqrt = TRUE)
  save(seurat__standard, 
       file = args[1], compress = FALSE)
} else {
  # save objects
  save(seurat__standard, 
       file = args[1], compress = FALSE)
}
