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
cell_info <- data.table::fread(rule$input$cell_info) # cell_info.tsv
cell_info$batch <- gsub(' ', '', cell_info$batch)
sample_meta <- data.table::fread(rule$input$sample_meta)

# set batch covariate for well data to NA, as any splits risks making the set too small
print('cell info import')
cell_info <- cell_info %>% 
  mutate(batch = case_when(study_accession == 'SRP125998' ~ paste0(study_accession, "_", Platform, '_NA'),
                           TRUE ~ batch)) %>% 
  mutate(batch = gsub(' ', '_', batch))
HVG_var_names =  scan(rule$input$HVG_var_names, what = 'character')
# cut down to approx max HVG length 
# assuming the list is in order of most var to least var
HVG_var_names = HVG_var_names[1:n_features] 
if (sum(is.na(HVG_var_names)) > 0) {
	HVG_var_names = HVG_var_names[!is.na(HVG_var_names)]
}
mito_geneids = scan(rule$input$mitogene_list, character(), sep= '\n') %>% 
  str_remove_all('\\.\\d+\\.$')
load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}

m <- load_rdata(rule$input$quant_file)



source(glue('{git_dir}/src/make_seurat_obj_functions.R') )

# load study exlucions list
exclude <- scan(glue('{git_dir}/data/exclusion.txt'), what = 'character') 
cell_info <- cell_info %>% filter(!study_accession %in% exclude, value %in% colnames(m))

# remove samples that are commented out in the srr_sample_file
keepers <- sample_meta %>% filter(!grepl('^#', sample_accession)) %>% pull(sample_accession)
cell_info <- cell_info %>% filter(sample_accession %in% keepers)

m <- m[, cell_info %>% pull(value)]
nFeature_RNA_cutoff = 500

if (set == 'universe'){
  seurat__standard <- make_seurat_obj(m[, cell_info %>% filter(!Source %in% c('Organoid', 'Cell Culture')) %>% pull(value)], cell_info, split.by = covariate, lengthCor = TRUE, dont_use_well_for_FVF = TRUE, mito_geneids = mito_geneids, nFeature_RNA_cutoff = 500)
} else if (set == 'universeHUMANHVG'){
  seurat__standard <- make_seurat_obj(m[, cell_info %>% filter(!Source %in% c('Organoid', 'Cell Culture')) %>% pull(value)], cell_info, split.by = covariate, lengthCor = TRUE, dont_use_well_for_FVF = TRUE, mito_geneids = mito_geneids, only_use_human_for_FVF = TRUE)
} else if (set == 'universeGIGAHVG'){
  seurat__standard <- make_seurat_obj(m[, cell_info %>% filter(!Source %in% c('Organoid', 'Cell Culture')) %>% pull(value)], cell_info, split.by = covariate, lengthCor = TRUE, dont_use_well_for_FVF = TRUE, mito_geneids = mito_geneids,  HVG = HVG_var_names)
} else if (set %in% c('mouse')){
  seurat__standard <- make_seurat_obj(m[, cell_info %>% filter(!Source %in% c('Organoid', 'Cell Culture')) %>% pull(value)], cell_info, split.by = covariate, lengthCor = TRUE, dont_use_well_for_FVF = TRUE, mito_geneids = mito_geneids)
} else if (set %in% c('chick', 'macaque', 'human')){
  seurat__standard <- make_seurat_obj(m[, cell_info %>% filter(!Source %in% c('Organoid', 'Cell Culture')) %>% pull(value)], cell_info, split.by = covariate, lengthCor = TRUE, dont_use_well_for_FVF = TRUE, mito_geneids = mito_geneids, keep_well = FALSE)
} else if (set == 'raw') {
  seurat__standard <- m
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
