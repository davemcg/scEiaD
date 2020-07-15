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
# custom combos
m_early <- m[,cell_info %>% filter(Age < 10) %>% pull(value)]
m_late <- m[,cell_info %>% filter(Age >= 10) %>% pull(value)]
m_test <- m[,sample(1:ncol(m), 10000)]
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


make_seurat_obj <- function(m, 
                            split.by = 'study_accession',
                            nfeatures = nfeatures,
							keep_well = TRUE,
							keep_droplet = TRUE,
							qumi = FALSE){
  well_m <- m[,cell_info %>% filter(value %in% colnames(m), !Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
  droplet_m <- m[,cell_info %>% filter(value %in% colnames(m), Platform %in% c('DropSeq', '10xv2', '10xv3')) %>% pull(value)]
  if (keep_well){
 	 seurat_well <- CreateSeuratObject(well_m)
  } 
  if (keep_droplet){
     seurat_droplet <- CreateSeuratObject(droplet_m)
  }
  
  # FILTER STEP!!!!
  # keep cells with < 10% mito genes, and more than 200 and less than 3000 detected genes for UMI
  # for well, drop the 3000 gene top end filter as there shouldn't be any droplets
  if (keep_well & !qumi){
	print('No QUMI')
    seurat_well <- subset(seurat_well, subset = nFeature_RNA > 200)
  } else if (keep_well && qumi) {
	print('QUMINORM!!')
    seurat_well <- subset(seurat_well, subset = nFeature_RNA > 200)
	qumi_counts <- quminorm(seurat_well@assays$RNA@counts)
  	seurat_well <- CreateSeuratObject(qumi_counts)
	seurat_well <- subset(seurat_well, subset = nFeature_RNA > 200)
	print('QUMI DONE!')
  }
  if (keep_droplet){
    seurat_droplet <- subset(seurat_droplet, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 )
  }
  # cells to keep
  if (!keep_well & keep_droplet){
    print('removing well')
	cells_to_keep <- row.names(seurat_droplet@meta.data)
  } else if (keep_well & !keep_droplet) {
    print('removing droplet')
	cells_to_keep <- row.names(seurat_well@meta.data)
  } else {
	print('keeping well and droplet')
    cells_to_keep <- c(row.names(seurat_droplet@meta.data), row.names(seurat_well@meta.data))
  }


  m_filter <- m[,cells_to_keep]
  
  seurat_m <- CreateSeuratObject(m_filter)
  seurat_m[["percent.mt"]] <- PercentageFeatureSet(seurat_m, pattern = "^MT-")
  seurat_m <- subset(seurat_m, subset = percent.mt < 10)
  
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
  seurat_m@meta.data$TechType <- left_join(seurat_m@meta.data %>% 
                                        row.names() %>% enframe(), 
                                      cell_info, by = 'value') %>% 
    pull(Platform)
  
  # scale data and regress
  seurat_m <- NormalizeData(seurat_m)
  # find var features

  seurat_m <- FindVariableFeatures(seurat_m, nfeatures = n_features, selection.method = 'vst')

  # don't use mito genes
  var_genes <- grep('^MT-', seurat_m@assays$RNA@var.features, value = TRUE, invert = TRUE)
  
  if (transform == 'standard'){
    print(paste0('Running lib.size and log correction, splitting by ', split.by))
    	
  	data <- seurat_m@assays$RNA@counts %>% t()
  	data <- data[,var_genes]
  	library_size <- Matrix::rowSums(data)
 	median_transcript_count <- stats::median(library_size)
  	data_norm <- median_transcript_count * data / library_size
  	data_norm <- t(data_norm)
  	data_norm <- log(data_norm + 1)
  	seurat_m@assays$RNA@scale.data <- data_norm %>% as.matrix()
     #seurat_m <- ScaleData(seurat_m,
      #                    features = var_genes,
       #                   split.by = split.by,
        #                  do.center = TRUE,
         #                 do.scale = TRUE,
          #                verbose = TRUE,
           #               vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
    
    seurat_m <- RunPCA(seurat_m, npcs = 100)
  }
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
  
  # remove sets with less than 500 cells, which will somehow(?) destory SCT - based integration performance
  low_n <- c()
  for (i in names(seurat_list)){
    if (ncol(seurat_list[[i]]) < 500){
      low_n <- c(low_n, i)
    }
  }
  if (length(low_n) > 0){
    seurat_list[low_n] <- NULL
  }
  
  study_data_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = n_features, verbose = FALSE)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = study_data_features, verbose = FALSE)
  
  # have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- RunPCA(x, features = study_data_features, verbose = FALSE)
  })
  list(seurat_list = seurat_list, study_data_features = study_data_features)
}

# scran normalization
scran_norm <- function(seurat_obj = seurat__standard, split.by = 'batch'){
  var_features <- grep('^MT-', seurat_obj@assays$RNA@var.features, invert =TRUE, value = TRUE)
  seurat_list <- SplitObject(seurat_obj, split.by = covariate)
  print('Beginning scran norm')
  # list of seurat objects
  sce_list <- list()
  for (obj in names(seurat_list)){
    print(obj)
    print(seurat_list[[obj]] %>% dim())
    if (ncol(seurat_list[[obj]]) > 50){
      sce_list[[obj]] <- SingleCellExperiment(assays = list(counts = as.matrix(x = seurat_list[[obj]]$RNA@data)))
      clusters = quickCluster(sce_list[[obj]], min.size=50)
      sce_list[[obj]] = computeSumFactors(sce_list[[obj]], cluster=clusters)
	  sce_list[[obj]] = scater::logNormCounts(sce_list[[obj]], log=TRUE)
      #sce_list[[obj]] = normalize(sce_list[[obj]], return_log = FALSE)
      print(summary(sizeFactors(sce_list[[obj]])))
      # seurat_list[[obj]]$RNA@data = as.sparse(log(x = assay(sce_list[[obj]], "normcounts") + 1))
	  seurat_list[[obj]]$RNA@data = assay(sce_list[[obj]], "logcounts") %>% as.sparse()
    } else {seurat_list[[obj]] <- NULL} # remove obj with less than 50 cells
  }
  # merge back into one seurat obj
  merged <- merge(x = seurat_list[[1]], y = seurat_list[2:length(x = seurat_list)])
  merged$RNA@scale.data = merged$RNA@data[var_features,] %>% as.matrix()
  merged@assays$RNA@var.features <- var_features
  # re do PCA
  merged <- RunPCA(merged, npcs = 100, features = var_features)
  merged
}

#' Performs L1 normalization on input data such that the sum of expression
#' values for each cell sums to 1, then returns normalized matrix to the metric
#' space using median UMI count per cell effectively scaling all cells as if
#' they were sampled evenly.
#' @param seurat_object
#' @return seurat_obj 
#' 2 dimensional array with normalized gene expression values
#' @import Matrix
#' @import dplyr
#'
#' @export
library.size.normalize <- function(seurat_obj, sqrt = FALSE, verbose=FALSE) {
  if (verbose) {
    message(paste0(
      "Normalizing library sizes for ",
      nrow(data), " cells"
    ))
  }
  vfeatures <- grep('^MT-', seurat_obj@assays$RNA@var.features, invert =TRUE, value = TRUE)
  data <- seurat_obj@assays$RNA@counts %>% t()
  data <- data[,vfeatures]
  library_size <- Matrix::rowSums(data)
  median_transcript_count <- stats::median(library_size)
  data_norm <- median_transcript_count * data / library_size
  data_norm <- t(data_norm)
  if (sqrt) {data_norm <- sqrt(data_norm)}
  seurat_obj@assays$RNA@scale.data <- data_norm %>% as.matrix()
  seurat_obj <- RunPCA(seurat_obj, features = vfeatures)
  seurat_obj
}

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

