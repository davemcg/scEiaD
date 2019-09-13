# combine datasets by species

# Use SCTransform to normalize the data, then RPCA to integrate
args <- commandArgs(trailingOnly = TRUE)

args <- c('quant/Mus_musculus/scTransformRPCA_anchor.seuratV3.Rdata','/home/mcgaugheyd/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv','references/gencode.vM22.metadata.MGI_tx_mapping.tsv','quant/Mus_musculus/counts.Rdata','quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')


library(Matrix)
library(tidyverse)
library(sva)
library(reticulate)
library(Seurat)
scanorama <- import('scanorama')
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)
downsample <- FALSE
stamp <- Sys.time() %>% gsub(' ', '__', .)

# load in metadata for study project merging, UMI correction, and gene name changing
metadata <- read_tsv(args[2])
tx <- read_tsv(args[3], col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

# well data
load(args[4])
# remove cells with > 10000 or < 1000
count <- count[,(diff(count@p) < 10000)]
count <- count[,(diff(count@p) > 1000)]

# sparse matrix files
rdata_files = args[5:length(args)]

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
  res_matrix <- res_matrix[,diff(res_matrix@p) < 6000]
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


# split by two groups
# early (< 10 days) and late (>10 days)
# Only one study (Clark ... Blackshaw is present for the early stage)
m_early <- m[,cell_info %>% filter(Age < 10) %>% pull(value)]
m_late <- m[,cell_info %>% filter(Age >= 10) %>% pull(value)]



make_seurat_obj <- function(m, normalize = FALSE, scale = FALSE){
  seurat_m <- CreateSeuratObject(m)
  seurat_m[["percent.mt"]] <- PercentageFeatureSet(seurat_m, pattern = "^MT-")
  # remove cells with > 10% mito genes
  seurat_m <- subset(seurat_m, subset = percent.mt < 10)
  # scale data and regress
  if (normalize){
    seurat_m <- NormalizeData(seurat_m)
  }
  # find var features
  seurat_m <- FindVariableFeatures(seurat_m, nfeatures = 3000, selection.method = 'vst')
  # var_genes <- grep('^MT-', seurat_m@assays$RNA@var.features, value = TRUE, invert = TRUE)
  if (scale){
    seurat_m <- ScaleData(seurat_m,
                          features = var_genes,
                          do.center = TRUE,
                          do.scale = TRUE,
                          vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
  }
  seurat_m@meta.data$batch <- cell_info %>% 
    filter(value %in% row.names(seurat_m@meta.data)) %>% 
    pull(batch)
  seurat_m@meta.data$Age <- cell_info %>% 
    filter(value %in% row.names(seurat_m@meta.data)) %>% 
    pull(Age)
  seurat_m
}


# split into two groups:
# < 10 days old (only one study)
# >= 10 days old (actually have independent studies)

seurat_m <- subset(seurat_m, subset = percent.mt < 10)

#s_data_list__early <- SplitObject(make_seurat_obj(m_early), split.by = 'batch')
s_data_list__late <- SplitObject(make_seurat_obj(m_late, normalize = TRUE), split.by = 'batch')

run_scanorama <- function(s_data_list, assay = 'RNA'){
  d <- list()
  g <- list()
  for (i in seq(1, length(s_data_list))){
    print(i);
    var_genes <- grep('^MT-', s_data_list[[i]]@assays$RNA@var.features, value = TRUE, invert = TRUE)
    d[[i]] <- t((s_data_list[[i]]@assays[[assay]]@data[var_genes,])) %>% as.matrix(); 
    d[[i]][is.na(d[[i]])] <- 0; 
    g[[i]] <- colnames(d[[i]]) 
  }
  integrated.corrected.data <- 
    scanorama$correct(d, g, return_dimred=FALSE, return_dense=TRUE)
  
  
  scan_cor <- Reduce(rbind, integrated.corrected.data[[1]])
  print('Running UMAP')
  umap <- uwot::umap(scan_cor, pca = 30, n_threads = 8)
  list(d = d, 
       g = g, 
       corrected_matrix = scan_cor, 
       umap_coordinates = umap)
}
#young <- run_scanorama(s_data_list__young)
old <- run_scanorama(s_data_list__late)
