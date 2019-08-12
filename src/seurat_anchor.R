# combine datasets by species

# Use SCTransform to normalize the data, then RPCA to integrate
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
plan(strategy = "multicore", workers = 12)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 80000 * 1024^2)
downsample <- FALSE

# load in metadata for study project merging, UMI correction, and gene name changing
metadata <- read_tsv(args[2])
tx <- read_tsv(args[3], col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

# well data
load(args[4])
# remove cells with > 10000 or < 1000
tpm <- tpm[,(diff(tpm@p) < 10000)]
tpm <- tpm[,(diff(tpm@p) > 1000)]

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

# grab available well samples
well_samples <- colnames(tpm)

# merge by project
study_sample <- metadata %>% 
  filter(sample_accession %in% c(droplet_samples, well_samples)) %>% 
  select(study_accession, Platform, sample_accession, Covariate) %>% 
  unique() %>% 
  mutate(study_accession = paste0(study_accession, '__', Platform, '__', Covariate)) %>%
  mutate(tech = case_when(sample_accession %in% droplet_samples ~ 'droplet',
                          TRUE ~ 'well')) %>% 
  arrange(study_accession)
# drop SRP161678 for now as it's giving weird errors right now
study_sample <- study_sample %>% filter(!grepl('SRP161678', study_accession))

# merge droplet and well data into one list of sparse matrices by study_sample$study_accession
# option to downsample to no more than 10,000 per study
study_data <- list()
for (i in unique(study_sample %>% pull(study_accession))){
  print(i)
  samples <- study_sample %>% filter(study_accession == i) %>% pull(sample_accession) 
  tech <- study_sample %>% filter(study_accession == i) %>% pull(tech) %>% head(1)
  if (tech == 'droplet'){
    study_data[[i]] <- do.call(cbind, sc_data[samples])
  }
  else {
    # remove samples that aren't in tpm (which means they failed upstream QC)
    samples <- samples[samples %in% well_samples]
    study_data[[i]] <- tpm[,samples]
    row.names(study_data[[i]]) <- row.names(study_data[[i]]) %>% toupper()
    row.names(study_data[[i]]) <- make.unique(row.names(study_data[[i]]))
    # remove na columns
    study_data[[i]] <- study_data[[i]][,colSums(is.na(study_data[[i]])) < 1]
  }
  if (downsample){
    sample_n <- min(ncol(study_data[[i]]),10000)
    cols <- sample(seq(1, sample_n))
    study_data[[i]] <- study_data[[i]][,cols]
  }
  # make seurat object, remove cells with few genes quantified
  study_data[[i]] <- CreateSeuratObject(study_data[[i]], project = i, min.cells = 250)
  # calc percentage mito genes
  study_data[[i]][["percent.mt"]] <- PercentageFeatureSet(study_data[[i]], pattern = "^MT-")
  study_data[[i]] <- SCTransform(study_data[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
  #study_data[[i]] <- RunPCA(study_data[[i]])
}

study_data_features <- SelectIntegrationFeatures(object.list = study_data, nfeatures = 3000, verbose = FALSE)
study_data <- PrepSCTIntegration(object.list = study_data, anchor.features = study_data_features, verbose = FALSE)

# have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
study_data <- lapply(X = study_data, FUN = function(x) {
  x <- RunPCA(x, features = study_data_features, verbose = FALSE)
})
save(study_data_features, study_data, file = 'study_data__emergency__2019_08_09.Rdata')

# try clear some memory
gc()
# nope, trying rpca now, running out of memory with the "CCT" reduction method
# 9,10 are "SRP075719__DropSeq__Batch1"     "SRP075719__DropSeq__Batch2"
anchors <- FindIntegrationAnchors(object.list = study_data, 
                                  normalization.method = 'SCT', 
                                  #reference = grep('SRP158081__10xv2', names(study_data)),
                                  scale = FALSE, 
                                  anchor.features = study_data_features, 
                                  reduction = "cca")

save(anchors, file = args[1])
