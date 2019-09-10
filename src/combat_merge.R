# combine datasets by species

# Use SCTransform to normalize the data, then RPCA to integrate
args <- commandArgs(trailingOnly = TRUE)

args <- c('quant/Mus_musculus/scTransformRPCA_anchor.seuratV3.Rdata','/home/mcgaugheyd/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv','references/gencode.vM22.metadata.MGI_tx_mapping.tsv','quant/Mus_musculus/counts.Rdata','quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')


library(Matrix)
library(tidyverse)
library(sva)
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

# create naive fully merged for combat sva correction
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

# find 2000 most variable genes
## skip mitochondrial
## for now use Seurat
library(Seurat)
seurat_m <- CreateSeuratObject(m)
seurat_m[["percent.mt"]] <- PercentageFeatureSet(seurat_m, pattern = "^MT-")
# remove cells with > 10% mito genes
seurat_m <- subset(seurat_m, subset = percent.mt < 10)
# find var features

var_genes <- grep('^MT-', seurat_m@assays$RNA@var.features, value = TRUE, invert = TRUE)

# scale data and regress
seurat_m <- NormalizeData(seurat_m)
seurat_m <- ScaleData(seurat_m, 
                      features = var_genes,
                      do.center = TRUE,
                      do.scale = TRUE,
                      vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seurat_m <- FindVariableFeatures(seurat_m, nfeatures = 5000)

small_m <- seurat_m@assays$RNA@scale.data[var_genes, ] %>% as.matrix()
small_m[is.na(small_m)] <- 0

# update cell_info to remove dropped cells
cell_info2 <- cell_info %>% filter(value %in% colnames(small_m)) %>% 
  mutate(age2 = case_when(Age < 10 ~ 0, Age >= 10 ~ 30, TRUE ~ Age),
         batch = paste(study_accession, Platform, Covariate, sep = '_'))
# combat
mod = model.matrix(object = ~as.factor(as.character(age2)),  
                   data = cell_info2)
combat_edata = ComBat(dat = small_m, 
                      batch=cell_info2$batch, 
                      mod = mod, 
                      par.prior=TRUE, 
                      prior.plots=FALSE)

# umap
umap <- uwot::umap(t(combat_edata), n_threads=8)