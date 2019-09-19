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
# load in metadata for study project merging, UMI correction, and gene name changing
transform = args[4]
metadata <- read_tsv(args[5])
tx <- read_tsv(args[6], col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

# well data
load(args[7])
# # remove cells with > 10000 or < 1000 detected genes
# count <- count[,(diff(count@p) < 10000)]
# count <- count[,(diff(count@p) > 1000)]

# extract species
species <- str_split(args[7], '/')[[1]][2]

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
cell_info <- cell_info %>% mutate(Age = case_when(Age > 100 ~ 30, TRUE ~ Age))
# save barcodes for labelling with published cell type assignment 
save(cell_info, file = paste0(species, '_cell_info.Rdata'))

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


# split into two groups:
# < 10 days old (only one study)
# >= 10 days old (actually have independent studies)


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

