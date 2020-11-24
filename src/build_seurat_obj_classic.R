# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger
args <- c('/data/swamyvs/scEiaD/rson_tmp/qjecvlvb.json', '/data/swamyvs/scEiaD/config.yaml')
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

gtf <- rtracklayer::readGFF('references/gtf/hs-homo_sapiens_anno.gtf.gz')
mito_geneids <- gtf %>% 
  filter(grepl('^MT-', ignore.case = T, gene_name)) %>% 
  pull(gene_id) %>% unique %>% str_remove_all('\\.\\d+$')

load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}

m <- load_rdata(rule$input$all_species_quant_file)

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


m_onlyDROPLET <- m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]

m_TABULA_DROPLET <- m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
m_onlyWELL <- m[,cell_info %>% filter(!Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]

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
  m <- m_early
  seurat__standard <- make_seurat_obj(m_early, split.by = covariate)
} else if (set == 'late'){
  print("Running Late")
  m <- m_late
  seurat__standard <- make_seurat_obj(m_late, split.by = covariate)
} else if (set == 'full'){
  print("Running Full")
  seurat__standard <- make_seurat_obj(m, split.by = covariate)
} else if (set == 'onlyDROPLET'){
  print("Running onlyDROPLET (remove well based)")
  m <- m_onlyDROPLET
  seurat__standard <- make_seurat_obj(m_onlyDROPLET, split.by = covariate, keep_well = FALSE)
}  else if (set == 'TabulaDroplet'){
  print("Running onlyDROPLET with Tabula Muris (no well)")
  m <- m_TABULA_DROPLET
  mito_geneids <- mito_geneids[mito_geneids%in% rownames(m_TABULA_DROPLET)]
  seurat__standard <- make_seurat_obj(m_TABULA_DROPLET, split.by = covariate, keep_well = FALSE,mito_geneids=mito_geneids)
} else if (set == 'onlyWELL' & transform == 'counts'){
  print("Running onlyWELL with quminorm (remove droplet based)") 
  m <- m_onlyWELL
  mito_geneids <- mito_geneids[mito_geneids%in% rownames(m_onlyWELL)]
  seurat__standard <- make_seurat_obj(m_onlyWELL, split.by = covariate, keep_droplet = FALSE, qumi = TRUE,mito_geneids=mito_geneids)
} else if (set == 'onlyWELL') {
  print("Running onlyWELL (remove droplet based)") 
  m <- m_onlyWELL
  mito_geneids <- mito_geneids[mito_geneids%in% rownames(m_onlyWELL)]
  seurat__standard <- make_seurat_obj(m_onlyWELL, split.by = covariate, keep_droplet = FALSE, mito_geneids=mito_geneids)
} else if (set == 'downsample'){
  print("Running downsample")
  m <- m_downsample
  seurat__standard <- make_seurat_obj(m_downsample, split.by = covariate)
} else if (set == 'universe'){
  seurat__standard <- make_seurat_obj(m, split.by = covariate, qumi = TRUE)
} else if (set %in% c('cones', 'hc', 'rgc', 'amacrine', 'mullerglia', 'bipolar', 'rods' )){
  m <- m_subset
  seurat__standard <- make_seurat_obj(m_subset, split.by = covariate, keep_well = FALSE)
} else if (set == 'raw') {
  seurat__standard <- m
} else if(set %in% c('Homo_sapiens', 'Mus_musculus')){ ##hacky way to loading in species specific data
  print(glue('loading {set} quant') )
  ## need to keep squaring off the data to keep uniformity for seurat
  
  if (set == 'Mus_musculus'){
    gtf <- rtracklayer::readGFF('references/gtf/mm-mus_musculus_anno.gtf.gz')
    mito_geneids <- gtf %>% 
      filter(grepl('^MT-', ignore.case = T, gene_name)) %>% 
      pull(gene_id) %>% unique %>% str_remove_all('\\.\\d+$')
    
  }
  
  species <- set 
  m_file <- str_replace_all(rule$input$all_species_quant_file, 'all_species_', glue('{species}/'))
  spec_m = load_rdata(m_file)
  seurat__standard =  make_seurat_obj(spec_m, split.by = covariate, keep_well = FALSE, mito_geneids=mito_geneids)
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
