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

# collapse by project
# merge by project
study_sample <- sra_metadata_extended %>% 
  filter(organism == 'Mus musculus', library_layout == 'PAIRED', UMI == 'YES') %>% 
  select(study_accession, sample_accession) %>% 
  unique() %>% 
  arrange(study_accession)

study_data <- list()
for (i in unique(study_sample$study_accession)){
  samples <- study_sample %>% filter(study_accession == i) %>% pull(sample_accession)
  study_data[[i]] <- do.call(cbind, sc_data[samples])
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


