# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger
args <- c('testing/test.yaml', '/data/swamyvs/scEiaD/config.yaml')
args <- commandArgs(trailingOnly = TRUE)
#Sys.setenv(SCIAD_CONFIG = '/data/swamyvs/scEiaD/config.yaml')

library(jsonlite)
library(yaml)
library(Matrix)
library(Matrix.utils)
library(tidyverse)
library(Seurat)
library(scran)
library(future)
library(quminorm)
library(glue)

rule <- read_json(args[1])
config <- read_yaml(args[2])


plan(strategy = "multicore", workers = 4)

git_dir=config$git_dir
working_dir=config$working_dir
setwd(working_dir)
print(working_dir)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)
print('load args')
set = rule$wildcards$partition # early, late, full, downsampled
covariate = rule$wildcards$covariate # study_accession, batch, etc.
transform = rule$wildcards$transform # SCT or standard seurat # mouse, mouse and macaque, mouse and macaque and human
n_features =as.numeric(rule$wildcards$n_features )
cell_info <-read_tsv(rule$input$cell_info) # cell_info.tsv
cell_info$batch <- gsub(' ', '', cell_info$batch)
# set batch covariate for well data to NA, as any splits risks making the set too small
print('cell info import')
cell_info <- cell_info %>% 
  mutate(batch = case_when(study_accession == 'SRP125998' ~ paste0(study_accession, "_", Platform, '_NA'),
                           TRUE ~ batch)) %>% 
  mutate(batch = gsub(' ', '_', batch))

mito_geneids = scan(rule$input$mitogene_list, character(), sep= '\n') %>% 
  str_remove_all('\\.\\d+\\.$')
load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}

m <- load_rdata(rule$input$all_species_quant_file)

print('Splitting time')
# # custom combos / sets

source(glue('{git_dir}/src/make_seurat_obj_functions.R') )

if (set == 'early'){
  print("Running Early")
  m_early <- m[,cell_info %>% filter(Age < 10) %>% pull(value)]
  seurat__standard <- make_seurat_obj(m_early,cell_info, split.by = covariate,mito_geneids=mito_geneids)
} else if (set == 'late'){
  print("Running Late")
  m_late <- m[,cell_info %>% filter(Age >= 10) %>% pull(value)]
  seurat__standard <- make_seurat_obj(m_late,cell_info, split.by = covariate,mito_geneids=mito_geneids)
} else if (set == 'full'){
  print("Running Full")
  seurat__standard <- make_seurat_obj(m,cell_info, split.by = covariate,mito_geneids=mito_geneids)
} else if (set == 'onlyDROPLET'){
  print("Running onlyDROPLET (remove well based)")
  m_onlyDROPLET <- m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]
  seurat__standard <- make_seurat_obj(m_onlyDROPLET, cell_info, split.by = covariate, keep_well = FALSE, mito_geneids = mito_geneids)
}  else if (set == 'TabulaDroplet'){
  print("Running onlyDROPLET with Tabula Muris (no well)")
  #m_TABULA_DROPLET <- m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
  m_TABULA_DROPLET <- m[,cell_info %>% filter(value %in% colnames(m), Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
  seurat__standard <- make_seurat_obj(m_TABULA_DROPLET,cell_info, split.by = covariate, keep_well = FALSE,mito_geneids=mito_geneids)
} else if (set == 'TabulaDropletLabelled'){
  print("Running onlyDROPLET labelled samples only with Tabula Muris (no well)")
  cells <- load_rdata('testing/TabulaDropletLabelled_barcodes.Rdata') %>% 
    .[.%in% colnames(m)]
  m <- m[,cells]
  seurat__standard <- make_seurat_obj(m,cell_info, split.by = covariate, keep_well = FALSE,mito_geneids=mito_geneids)
}else if (set == 'onlyWELL' & transform == 'counts'){
  print("Running onlyWELL with quminorm (remove droplet based)") 
  m_onlyWELL <- m[,cell_info %>% filter(!Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]
  seurat__standard <- make_seurat_obj(m_onlyWELL,cell_info, split.by = covariate, keep_droplet = FALSE, qumi = TRUE,mito_geneids=mito_geneids)
} else if (set == 'onlyWELL') {
  print("Running onlyWELL (remove droplet based)") 
  m_onlyWELL <- m[,cell_info %>% filter(!Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]
  seurat__standard <- make_seurat_obj(m_onlyWELL,cell_info, split.by = covariate, keep_droplet = FALSE, mito_geneids=mito_geneids)
} else if (set == 'downsample'){
  print("Running downsample")
  downsample_samples <- 
    cell_info %>% 
    group_by(batch) %>% 
    sample_n(2000, replace = TRUE) %>% 
    unique() %>% 
    pull(value)
  m_downsample <- m[,downsample_samples]
  seurat__standard <- make_seurat_obj(m_downsample,cell_info, split.by = covariate,mito_geneids=mito_geneids)
} else if (set == 'universe'){
  seurat__standard <- make_seurat_obj(m, cell_info, split.by = covariate, lengthCor = TRUE, dont_use_well_for_FVF = TRUE, mito_geneids = mito_geneids)
} else if (set == 'universeX'){
  seurat__standard <- make_seurat_obj(m, cell_info, split.by = covariate, lengthCor = TRUE, only_use_human_for_FVF = TRUE, mito_geneids = mito_geneids)
} else if (set %in% c('cones', 'hc', 'rgc', 'amacrine', 'mullerglia', 'bipolar', 'rods' )){
  # no circular dependencies!
  precursors <- c('AC/HC_Precurs', 'Early RPCs', 'Late RPCs', 'Neurogenic Cells', 'Photoreceptor Precursors', 'RPCs')
  if (set == 'cones') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Cones', precursors)) %>% pull(value)]
  }
  if (set == 'bipolar') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Bipolar Cells', precursors)) %>% pull(value)]
  }
  if (set == 'rods') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Rods', precursors)) %>% pull(value)]
  }
  if (set == 'mullerglia') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Muller Glia', precursors)) %>% pull(value)]
  }
  if (set == 'amacrine') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Amacrine Cells', precursors)) %>% pull(value)]
  }
  if (set == 'rgc') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Retinal Ganglion Cells', precursors)) %>% pull(value)]
  }
  if (set == 'hc') {
  	m_subset = m[, cell_info %>% filter(CellType %in% c('Horizontal Cells', precursors)) %>% pull(value)]
  }
  seurat__standard <- make_seurat_obj(m_subset, cell_info, split.by = covariate, keep_well = FALSE, mito_geneids = mito_geneids)
} else if (set == 'raw') {
  seurat__standard <- m
} else if(set %in% c('Homo_sapiens', 'Mus_musculus')){ ##hacky way to loading in species specific data
  print(glue('loading {set} quant') )
  ## need to keep squaring off the data to keep uniformity for seurat
  
  if (set == 'Mus_musculus'){
    mito_geneids=scan('references/mito_genes/mm-mus_musculus_mitogenes.txt', character(), sep = '\n') %>% 
      str_remove_all('\\.\\d+\\.$')
    
  }
  
  species <- set 
  m_file <- str_replace_all(rule$input$all_species_quant_file, 'all_species_', glue('{species}/'))
  spec_m = load_rdata(m_file)
  seurat__standard =  make_seurat_obj(spec_m,cell_info, split.by = covariate, keep_well = FALSE, mito_geneids=mito_geneids)
} 



if (transform == 'SCT'){
  s_data_list<- SplitObject(seurat__standard, split.by = covariate)
  seurat__SCT <- seurat_sct(s_data_list)
  save(seurat__SCT, 
       file = rule$output$seurat, compress = FALSE)
} else if (transform == 'scran'){
  seurat__standard <- scran_norm(seurat__standard, split.by = covariate )
  save(seurat__standard, 
       file = rule$output$seurat, compress = FALSE)
} else if (transform == 'libSize'){
  seurat__standard <- library.size.normalize(seurat__standard)
  save(seurat__standard, 
       file = rule$output$seurat, compress = FALSE)
} else if (transform == 'sqrt') {
  seurat__standard <- library.size.normalize(seurat__standard, sqrt = TRUE)
  save(seurat__standard, 
       file = rule$output$seurat, compress = FALSE)
} else {
  # save objects
  save(seurat__standard, 
       file = rule$output$seurat, compress = FALSE)
}
