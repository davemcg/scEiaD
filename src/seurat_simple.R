library(Seurat)
library(tidyverse)
library(Matrix)
library(jsonlite)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
rule <- read_json(args[1])
config <- read_yaml(args[2])

print(rule)
print(config)

git_dir=config$git_dir
working_dir=config$working_dir
setwd(working_dir)
print(working_dir)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 500000 * 1024^2)
sample = rule$wildcard$sample #e.g. SRX8884952
n_features =as.numeric(rule$wildcards$n_features )
RNA_count_min = 450
res = 0.8
n_dims = as.numeric(rule$wildcards$dims)
cell_info <- data.table::fread(rule$input$cell_info) # cell_info.tsv
cell_info$batch <- gsub(' ', '', cell_info$batch)
load(rule$input$labelled_cells)
sceiad_meta <- fst::read_fst(rule$input$scEiaD_meta)
if (rule$wildcards$norm == 'SCT'){
	SCT <- TRUE
} else {
	SCT <- FALSE
}
print('Load Mito')
mito_files <- c('references/mito_genes/gg-gallus_gallus_mitogenes.txt',
					'references/mito_genes/mm-mus_musculus_mitogenes.txt',
					'references/mito_genes/hs-homo_sapiens_mitogenes.txt')
mito_l <- list()
for (i in mito_files){
	mito_l[[i]] <- scan(i, character(), sep= '\n') %>%
  str_remove_all('\\.\\d+\\.$')
}
mito_geneids = mito_l %>% unlist() %>% unique() 

load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}
print('Load Big Matrix')
m <- load_rdata(rule$input$count_matrix)

# update batch
srl <-  read_tsv('~/git/scEiaD/data/sample_run_layout_organism_tech_biosample_organ_2022_02_18.tsv')
cell_info <- cell_info %>% select(-Covariate) %>% left_join(srl %>% select(sample_accession, Covariate) %>% unique()) 
cell_info <- cell_info %>% mutate(batch = paste0(study_accession, '_', Platform, '_', Covariate))

# cut down m to just the batch given
cells_to_retain <- cell_info %>% filter(batch == sample) %>% pull(value)

m_sample <- m[,cells_to_retain]

print("Seurat time")
seurat <- CreateSeuratObject(m_sample)
seurat <- subset(seurat, subset = nFeature_RNA > RNA_count_min)
mito_geneids_present <- mito_geneids[mito_geneids %in% rownames(seurat)]
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mito_geneids_present)

seurat@meta.data$batch <- left_join(seurat@meta.data %>%
                                          row.names() %>% enframe(),
                                        cell_info, by = 'value') %>%
    pull(batch)
seurat@meta.data$study_accession <- left_join(seurat@meta.data %>%
                                                    row.names() %>% enframe(),
                                                  cell_info, by = 'value') %>%
    pull(study_accession)
seurat@meta.data$Age <- left_join(seurat@meta.data %>%
                                        row.names() %>% enframe(),
                                      cell_info, by = 'value') %>%
    pull(Age)
seurat@meta.data$TechType <- left_join(seurat@meta.data %>%
                                        row.names() %>% enframe(),
                                      cell_info, by = 'value') %>%
    pull(Platform)
seurat@meta.data$CellType <- left_join(seurat@meta.data %>%
                                        row.names() %>% enframe(),
                                      cell_info_labels, by = 'value') %>%
    pull(CellType)
# MT filtering
# 10xv3 has approx 2x the mito counts of 10xv2, so giving that platform a diff cutoff
# https://kb.10xgenomics.com/hc/en-us/articles/360026501692-Do-we-see-a-difference-in-expression-profile-of-3-Single-Cell-v3-chemistry-as-compared-to-v2-chemistry
mt_cells_keep <- (seurat@meta.data$percent.mt < 10 & seurat@meta.data$TechType != '10xv3') | (seurat@meta.data$percent.mt < 20 & seurat@meta.data$TechType == '10xv3')
seurat <- seurat[, mt_cells_keep]
if (SCT){
	print("SCT NORM")
	seurat <- SCTransform(seurat, variable.features.n = n_features + 1000, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE, method = "glmGamPoi")
	n_dims = n_dims + 10
 } else { 
	print("STANDARD NORM")		
	seurat <- NormalizeData(seurat)
	seurat <- FindVariableFeatures(seurat, nfeatures = n_features, selection.method = 'vst')
	seurat <- ScaleData(seurat)
}
seurat <- RunPCA(seurat, npcs = 50)
seurat <- FindNeighbors(seurat, dims = 1:n_dims)
seurat <- FindClusters(seurat, resolution = res)
seurat <- RunUMAP(seurat, dims = 1:n_dims)

# add scEiaD meta
sample_meta <- seurat@meta.data
sample_meta <- sample_meta %>% as_tibble(rownames = 'Barcode')  %>% left_join(sceiad_meta %>% select(one_of(colnames(sceiad_meta)[!colnames(sceiad_meta) %in% colnames(sample_meta)]))) %>% data.frame()
row.names(sample_meta) <- sample_meta$Barcode
seurat@meta.data <- sample_meta

save(seurat, file = rule$output$seurat)

