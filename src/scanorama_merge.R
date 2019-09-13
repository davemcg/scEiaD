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
  var_genes <- grep('^MT-', seurat_m@assays$RNA@var.features, value = TRUE, invert = TRUE)
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

# seurat_early <- make_seurat_obj(m_early, normalize = TRUE, scale = TRUE)
# s_data_list__early <- SplitObject(seurat_early, split.by = 'batch')
seurat_late <- make_seurat_obj(m_late, normalize = TRUE, scale = TRUE)
s_data_list__late <- SplitObject(seurat_late, split.by = 'batch')


run_scanorama <- function(s_data_list, assay = 'RNA', num_dims = 30, run_umap = FALSE){
  # scanorama can return "integrated" and/or "corrected" data
  # authors say that the "integrated" data is a low-dimension (100) representation
  # of the integration, which is INTENDED FOR PCA/tSNE/UMAP!!!
  # the corrected data returns all of the sample x gene matrix with batch
  # corrected values
  d <- list()
  g <- list()
  for (i in seq(1, length(s_data_list))){
    print(i);
    var_genes <- grep('^MT-', s_data_list[[i]]@assays$RNA@var.features, value = TRUE, invert = TRUE)
    d[[i]] <- t((s_data_list[[i]]@assays[[assay]]@scale.data[var_genes,])) %>% as.matrix(); 
    d[[i]][is.na(d[[i]])] <- 0; 
    g[[i]] <- colnames(d[[i]]) 
  }

  integrated.corrected.data <- 
    scanorama$correct(d, g , return_dimred=TRUE, return_dense=TRUE)
  # relabel values as reticulate/scanorama don't return the row or column names
  for (i in seq(1, length(d))){
    # first is in the integrated data (dim reduced for UMAP, etc)
    row.names(integrated.corrected.data[[1]][[i]]) <- row.names(d[[i]])
    # second is the corrected gene expression values
    row.names(integrated.corrected.data[[2]][[i]]) <- row.names(d[[i]])
    colnames(integrated.corrected.data[[2]][[i]]) <- colnames(d[[i]])
  }

  scanorama_int_matrix <- Reduce(rbind, integrated.corrected.data[[1]])
  if (run_umap){
  print('Running UMAP')
  umap <- uwot::umap(scanorama_int_matrix[,1:num_dims], n_threads = 8)
  names <- d %>% map(row.names) %>% flatten_chr()
  row.names(umap) <- names
  colnames(umap) <- c('UMAP_1', 'UMAP_2')
  } else {umap = 'not run'}
  list(d = d, 
       g = g, 
       scanorama = integrated.corrected.data,
       scanorama_integrated_matrix = scanorama_int_matrix, 
       umap_coordinates = umap)
  
}
#early <- run_scanorama(s_data_list__young)
late <- run_scanorama(s_data_list__late)


# put scanorama low dim reduction into seurat obj
scanorama_mnn <- late$scanorama_integrated_matrix
colnames(scanorama_mnn) <- paste0("scanorama_", 1:ncol(scanorama_mnn))
seurat_late[["scanorama"]] <- CreateDimReducObject(embeddings = scanorama_mnn, key = "scanorama_", assay = DefaultAssay(seurat_late))

seurat_late <- RunUMAP(seurat_late, dims = 1:20, reduction = 'scanorama', reduction.key = 'scanoramaUMAP_')





# bbknn

pca <- irlba(seurat_late@assays$RNA@scale.data %>% as.matrix(), nv = 50)
pca <- pca$v
batch <- seurat_late@meta.data$batch %>% as.character() %>% as.factor()

adata = anndata$AnnData(X=pca, obs=batch)
sc$tl$pca(adata)
adata$obsm$X_pca = pca
bbknn$bbknn(adata,batch_key=0)
sc$tl$umap(adata)
umap = py_to_r(adata$obsm$X_umap)

