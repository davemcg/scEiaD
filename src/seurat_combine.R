# combine datasets by species

# first with 10X/UMI and mouse
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
plan(strategy = "multicore", workers = 16)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 40000 * 1024^2)

# load well (not UMI/droplet) data
load(args[2])

rdata_files = args[3:length(args)]
sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  sc_data[[sample_accession]] <- seurat_obj
}

sc_data_features <- SelectIntegrationFeatures(object.list = sc_data, nfeatures = 3000, verbose = FALSE)
sc_data <- PrepSCTIntegration(object.list = sc_data, anchor.features = sc_data_features, verbose = FALSE)

# have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
sc_data <- lapply(X = sc_data, FUN = function(x) {
  x <- RunPCA(x, features = sc_data_features, verbose = FALSE)
})

# nope, trying rpca now, running out of memory with the "CCT" reduction method
anchors <- FindIntegrationAnchors(object.list = sc_data, 
                                  normalization.method = 'SCT', 
                                  scale = FALSE, 
                                  anchor.features = sc_data_features, 
                                  reduction = "rpca")

sc_data_integrated <- IntegrateData(anchorset = anchors, normalization.method = 'SCT', verbose = TRUE)

# PCA, UMAP
sc_data_integrated <- RunPCA(sc_data_integrated, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
sc_data_integrated <- RunUMAP(sc_data_integrated, dims = 1:75, min.dist = 0.75)

pdf('merged.pdf')
print(DimPlot(sc_data_integrated, group.by = 'orig.ident'))
dev.off()

save(sc_data_integrated, file = args[1])

















## study level
tx <- read_tsv('references/gencode.vM22.metadata.MGI_tx_mapping.tsv', col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')
rdata_files <- c('quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')
sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  row.names(res_matrix) <- row.names(res_matrix) %>% 
    enframe(value = 'id') %>% 
    left_join(., tx, by = 'id') %>% 
    pull(gene)
  colnames(res_matrix) <- make.unique(colnames(res_matrix))
  colnames(res_matrix) <- paste0(colnames(res_matrix), "_", sample_accession)
  row.names(res_matrix) <- make.unique(row.names(res_matrix))
  sc_data[[sample_accession]] <- res_matrix
}

# merge by project
study_sample <- sra_metadata_extended %>% filter(organism == 'Mus musculus', library_layout == 'PAIRED', UMI == 'YES') %>% select(study_accession, biosample_attribute_recs, sample_accession) %>% unique() %>% arrange(study_accession)
study_sample <- study_sample %>% mutate(study_accession = case_when(grepl('Mixed', biosample_attribute_recs) ~ paste0(study_accession, '_P14_Mixed_CD1_C57Bl6'),
                                                                    TRUE ~ study_accession))
study_sample 
study_data <- list()
for (i in unique(study_sample$study_accession)){
  samples <- study_sample %>% filter(study_accession == i) %>% pull(sample_accession)
  study_data[[i]] <- do.call(cbind, sc_data[samples])
  study_data[[i]] <- CreateSeuratObject(study_data[[i]])
  study_data[[i]] <- SCTransform(study_data[[i]])
  study_data[[i]] <- RunPCA(study_data[[i]])
}

study_data_features <- SelectIntegrationFeatures(object.list = study_data, nfeatures = 3000, verbose = FALSE)
study_data <- PrepSCTIntegration(object.list = study_data, anchor.features = study_data_features, verbose = FALSE)

# have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
study_data <- lapply(X = study_data, FUN = function(x) {
  x <- RunPCA(x, features = study_data_features, verbose = FALSE)
})

# nope, trying rpca now, running out of memory with the "CCT" reduction method
anchors <- FindIntegrationAnchors(object.list = study_data, 
                                  normalization.method = 'SCT', 
                                  scale = FALSE, 
                                  anchor.features = study_data_features, 
                                  reduction = "rpca")

study_data_integrated <- IntegrateData(anchorset = anchors, normalization.method = 'SCT', verbose = TRUE)

# PCA, UMAP
study_data_integrated <- RunPCA(study_data_integrated, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
study_data_integrated <- RunUMAP(study_data_integrated, dims = 1:75, min.dist = 0.75)

pdf('merged.pdf')
print(DimPlot(study_data_integrated, group.by = 'orig.ident'))
dev.off()

save(study_data_integrated, file = args[1])

