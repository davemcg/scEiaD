# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)

#args <- c('seurat_obj/Mus_musculus__standard_and_SCT__late__batch.seuratV3.Rdata','late','batch','/home/mcgaugheyd/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv','references/gencode.vM22.metadata.MGI_tx_mapping.tsv','quant/Mus_musculus/counts.Rdata','quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

library(Matrix)
library(tidyverse)
library(Seurat)
library(scran)
library(future)
library(quminorm)
plan(strategy = "multicore", workers = 4)
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

nfeatures = args[7] %>% as.numeric()
rdata_files = args[8:length(args)]
rdata <- list()
for (i in rdata_files){
  load(i)
  rdata[[i]] <- m
}
if (combination == 'Mus_musculus'){
  file <- grep('Mus_mus', names(rdata), value = TRUE)
  m <- rdata[[file]]
} else if (combination == 'MacaMusHomoMuris') {
  #  mito naming in macaque missing the 'MT-' part....
  # also call CO1 and CO2 -> COX1 and COX2
  # sigh
  mf_gene <- rdata[["quant/Macaca_fascicularis/full_sparse_matrix.Rdata"]] %>% row.names()
  mf_gene <- mf_gene %>% enframe() %>% mutate(value = case_when(grepl('^ND\\d', value) ~ paste0('MT-', value),
                                                     value == 'COX1' ~ 'MT-CO1',
                                                     value == 'COX2' ~ 'MT-CO2',
                                                     value == 'COX3' ~ 'MT-CO3',
                                                     value == 'ATP6' ~ 'MT-ATP6',
                                                     value == 'ATP8' ~ 'MT-ATP8',
                                                     TRUE ~ value ))
  row.names(rdata[["quant/Macaca_fascicularis/full_sparse_matrix.Rdata"]]) <- mf_gene$value
  load('tabula_muris_combined.Rdata') 
  # update gene symbols
  mgi <- read_tsv('~/git/massive_integrated_eye_scRNA/data/MGIBatchReport_20200701_114356.txt')
  new_names <- row.names(facs_mat) %>% toupper() %>% as_tibble() %>% left_join(mgi, by = c('value' = 'Input')) %>% mutate(nname = case_when(is.na(Symbol) ~ value, TRUE ~ toupper(Symbol))) %>% group_by(value) %>% summarise(nname = head(nname, 1))
  row.names(facs_mat) <- new_names$nname
  new_names <- row.names(droplet_mat) %>% toupper() %>% as_tibble() %>% left_join(mgi, by = c('value' = 'Input')) %>% mutate(nname = case_when(is.na(Symbol) ~ value, TRUE ~ toupper(Symbol))) %>% group_by(value) %>% summarise(nname = head(nname, 1))
  row.names(droplet_mat) <- new_names$nname

  shared_genes <- rdata %>% map(row.names) %>% purrr::reduce(intersect)
  missing_in_muris <- shared_genes[!(shared_genes %in% row.names(facs_mat))]

  empty <- Matrix(0, nrow = length(missing_in_muris), ncol = ncol(facs_mat), sparse = TRUE)
  row.names(empty) <- missing_in_muris
  facs_mat <- rbind(facs_mat, empty)

  empty <- Matrix(0, nrow = length(missing_in_muris), ncol = ncol(droplet_mat), sparse = TRUE)
  row.names(empty) <- missing_in_muris
  droplet_mat <- rbind(droplet_mat, empty)

  if (set == 'onlyDROPLET'){
	print('Adding Tabula Muris Droplet')
    rdata[['droplet_muris']] <- droplet_mat
    cell_info <- bind_rows(cell_info, 
							droplet_meta %>% mutate( 
								batch = paste0('TabulaMuris_', `mouse.id`),
								organism = 'Mus musculus',
 								Tissue = tissue,
								sample_accession = paste0('TabulaMuris_', tissue),
								study_accession = 'TabulaMuris',
								UMI = 'YES',
								library_layout = 'PAIRED',
								integration_group = 'Late',
								Platform = '10xv2') %>%
							select(value, batch:Platform))
  } else if (set == 'onlyWELL') { 
	print('Adding Tabula Muris FACS')
    rdata[['facs_muris']] <- facs_mat
    cell_info <- bind_rows(cell_info, 
							facs_meta %>% mutate( 
								batch = paste0('TabulaMuris_', `mouse.id`),
								organism = 'Mus musculus',
 								Tissue = tissue,
								sample_accession = paste0('TabulaMuris_', tissue),
								study_accession = 'TabulaMuris',
								UMI = 'NO',
								library_layout = 'PAIRED',
								integration_group = 'Late',
								Platform = 'SMARTSeq_v2') %>%
							select(value, batch:Platform))
  }
  file_cut_down <- list()
  mito_list <- list()
  for (i in names(rdata)){
    file_cut_down[[i]] <- rdata[[i]][shared_genes,]
  }
  m <- file_cut_down %>% purrr::reduce(cbind)
} else {
 
  #  mito naming in macaque missing the 'MT-' part....
  # also call CO1 and CO2 -> COX1 and COX2
  # sigh
  mf_gene <- rdata[["quant/Macaca_fascicularis/full_sparse_matrix.Rdata"]] %>% row.names()
  mf_gene <- mf_gene %>% enframe() %>% mutate(value = case_when(grepl('^ND\\d', value) ~ paste0('MT-', value),
                                                     value == 'COX1' ~ 'MT-CO1',
                                                     value == 'COX2' ~ 'MT-CO2',
                                                     value == 'COX3' ~ 'MT-CO3',
                                                     value == 'ATP6' ~ 'MT-ATP6',
                                                     value == 'ATP8' ~ 'MT-ATP8',
                                                     TRUE ~ value ))
  row.names(rdata[["quant/Macaca_fascicularis/full_sparse_matrix.Rdata"]]) <- mf_gene$value

  shared_genes <- rdata %>% map(row.names) %>% purrr::reduce(intersect)
  file_cut_down <- list()
  mito_list <- list()
  for (i in names(rdata)){
    file_cut_down[[i]] <- rdata[[i]][shared_genes,]
  }
  m <- file_cut_down %>% purrr::reduce(cbind)
}

print('Splitting time')
# custom combos / sets
m_early <- m[,cell_info %>% filter(Age < 10) %>% pull(value)]
m_late <- m[,cell_info %>% filter(Age >= 10) %>% pull(value)]
m_test <- m[,sample(1:ncol(m), 10000)]


precursors <- c('AC/HC_Precurs', 'Early RPCs', 'Late RPCs', 'Neurogenic Cells', 'Photoreceptor Precursors', 'RPCs')
load('/data/OGVFB_BG/scEiaD/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umapFilter.predictions.Rdata')
if (set == 'cones') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Cones', precursors)) %>% pull(Barcode)]
}
if (set == 'bipolar') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Bipolar Cells', precursors)) %>% pull(Barcode)]
}
if (set == 'rods') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Rods', precursors)) %>% pull(Barcode)]
}
if (set == 'mullerglia') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Muller Glia', precursors)) %>% pull(Barcode)]
}
if (set == 'amacrine') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Amacrine Cells', precursors)) %>% pull(Barcode)]
}
if (set == 'rgc') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Retinal Ganglion Cells', precursors)) %>% pull(Barcode)]
}
if (set == 'hc') {
	m_subset = m[, umap %>% filter(CellType_predict %in% c('Horizontal Cells', precursors)) %>% pull(Barcode)]
}


m_onlyDROPLET <-  m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]
m_TABULA_DROPLET <- m[,cell_info %>% filter(Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
m_onlyWELL <- m[,cell_info %>% filter(!Platform %in% c('DropSeq', '10xv2', '10xv3'), study_accession != 'SRP131661') %>% pull(value)]

downsample_samples <- 
  cell_info %>% 
  group_by(batch) %>% 
  sample_n(2000, replace = TRUE) %>% 
  unique() %>% 
  pull(value)
m_downsample <- m[,downsample_samples]

source('/home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/make_seurat_obj_functions.R')

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
} else if (set == 'onlyWELL' & transform == 'counts'){
  print("Running onlyWELL with quminorm (remove droplet based)") 
  seurat__standard <- make_seurat_obj(m_onlyWELL, split.by = covariate, keep_droplet = FALSE, qumi = TRUE)
} else if (set == 'onlyWELL') {
  print("Running onlyWELL (remove droplet based)") 
  seurat__standard <- make_seurat_obj(m_onlyWELL, split.by = covariate, keep_droplet = FALSE)
} else if (set == 'downsample'){
  print("Running downsample")
  seurat__standard <- make_seurat_obj(m_downsample, split.by = covariate)
} else if (set == 'universe'){
  seurat__standard <- make_seurat_obj(m, split.by = covariate, qumi = TRUE)
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

