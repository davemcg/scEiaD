# build seurat obj with the "classic" findvariablegenes -> normalize -> scaledata processing
# output use for integration with various algorithms
# scanorama, CCT, harmony, liger

args <- commandArgs(trailingOnly = TRUE)

#args <- c('seurat_obj/Mus_musculus__standard_and_SCT__late__batch.seuratV3.Rdata','late','batch','/home/mcgaugheyd/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv','references/gencode.vM22.metadata.MGI_tx_mapping.tsv','quant/Mus_musculus/counts.Rdata','quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

library(Matrix)
library(tidyverse)
library(Seurat)
library(future)
plan(strategy = "multicore", workers = 4)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)

set = args[2] # early, late, full
covariate = args[3] # study_accession, batch, etc.
transform = args[4] # SCT or standard seurat
combination = args[5] # mouse, mouse and macaque, mouse and macaque and human
cell_info <- read_tsv(args[6]) # cell_info.tsv
rdata_files = args[7:length(args)]
rdata <- list()
for (i in rdata_files){
  load(i)
  rdata[[i]] <- m
}
if (combination == 'Mus_musculus'){
  file <- grep('Mus_mus', names(rdata), value = TRUE)
  m <- rdata[[file]]
} else if (combination == 'Mus_musculus_Macaca_fascicularis'){
  file  <- grep('Mus|Maca', names(rdata), value = TRUE)
  shared_genes <- rdata[file] %>% map(row.names) %>% reduce(intersect)
  file_cut_down <- list()
  for (i in file){
    file_cut_down[[i]] <- rdata[[i]][shared_genes,]
  }
  m <- file_cut_down %>% reduce(cbind)
} else {
  shared_genes <- rdata %>% map(row.names) %>% reduce(intersect)
  file_cut_down <- list()
  for (i in names(rdata)){
    file_cut_down[[i]] <- rdata[[i]][shared_genes,]
  }
  m <- file_cut_down %>% reduce(cbind)
}


# split by two groups
# early (< 10 days) and late (>10 days)
# Only one study (Clark ... Blackshaw is present for the early stage)
m_early <- m[,cell_info %>% filter(Age < 10) %>% pull(value)]
m_late <- m[,cell_info %>% filter(Age >= 10) %>% pull(value)]
m_test <- m[,sample(1:ncol(m), 10000)]


make_seurat_obj <- function(m, 
                            split.by = 'study_accession'){
  seurat_m <- CreateSeuratObject(m)
  seurat_m[["percent.mt"]] <- PercentageFeatureSet(seurat_m, pattern = "^MT-")
  # FILTER STEP!!!!
  # keep cells with < 10% mito genes, and more than 200 and less than 3000 detected genes
  seurat_m <- subset(seurat_m, subset = percent.mt < 10 & nFeature_RNA > 200 & nFeature_RNA < 3000 )
  seurat_m@meta.data$batch <- left_join(seurat_m@meta.data %>% 
                                          row.names() %>% enframe(), 
                                        cell_info, by = 'value') %>% 
    pull(batch)
  seurat_m@meta.data$study_accession <- left_join(seurat_m@meta.data %>% 
                                                    row.names() %>% enframe(), 
                                                  cell_info, by = 'value') %>% 
    pull(study_accession)
  seurat_m@meta.data$Age <- left_join(seurat_m@meta.data %>% 
                                        row.names() %>% enframe(), 
                                      cell_info, by = 'value') %>% 
    pull(Age)
  # scale data and regress
  seurat_m <- NormalizeData(seurat_m)
  # find var features
  seurat_m <- FindVariableFeatures(seurat_m, nfeatures = 2000, selection.method = 'vst')
  # don't use mito genes
  var_genes <- grep('^MT-', seurat_m@assays$RNA@var.features, value = TRUE, invert = TRUE)
  
  print(paste0('Running scale, splitting by ', split.by))
  seurat_m <- ScaleData(seurat_m,
                        features = var_genes,
                        split.by = split.by,
                        do.center = TRUE,
                        do.scale = TRUE,
                        verbose = TRUE,
                        vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
  
  seurat_m <- RunPCA(seurat_m, npcs = 100)
  seurat_m
}


# build SCT based seurat obj
seurat_sct <- function(seurat_list){
  ## tryCatch for SCT
  trySCTransform <- function(x){
    tryCatch(
      expr = {
        message("Successful SCTransform")
        return(SCTransform(x, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")))
      },
      error = function(e){
        message("Failed SCTransform")
        print(e)
        return(NULL)
      }
    )
  }
  for (i in names(seurat_list)){
    DefaultAssay(seurat_list[[i]]) <- 'RNA'
    seurat_list[[i]] <- trySCTransform(seurat_list[[i]])
  }
  
  # remove sets with less than 200 cells, which will fail integration
  low_n <- c()
  for (i in names(seurat_list)){
    if (ncol(seurat_list[[i]]) < 200){
      low_n <- c(low_n, i)
    }
  }
  if (length(low_n) > 0){
    seurat_list[low_n] <- NULL
  }
  
  study_data_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000, verbose = FALSE)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = study_data_features, verbose = FALSE)
  
  # have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- RunPCA(x, features = study_data_features, verbose = FALSE)
  })
  
  list(seurat_list = seurat_list, study_data_features = study_data_features)
  
}

if (set == 'early'){
  print("Running Early")
  seurat__standard <- make_seurat_obj(m_early, split.by = covariate)
} else if (set == 'late'){
  print("Running Late")
  seurat__standard <- make_seurat_obj(m_late, split.by = covariate)
} else {
  print("Running Full")
  seurat__standard <- make_seurat_obj(m, split.by = covariate)
}

if (transform == 'SCT'){
  s_data_list<- SplitObject(seurat__standard, split.by = covariate)
  seurat__SCT <- seurat_sct(s_data_list)
  save(seurat__SCT, 
       file = args[1], compress = FALSE)
} else {
  # save objects
  save(seurat__standard, 
       file = args[1], compress = FALSE)
}

