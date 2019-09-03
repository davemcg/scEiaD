# combine datasets by species

# Use SCTransform to normalize the data, then RPCA to integrate
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
plan(strategy = "multicore", workers = 8)
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
count <- tpm
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
  samples <- c(droplet_samples, sample_accession)
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
well_samples <- colnames(count)

# merge by project, platform, covariate, and integration_group
study_sample <- metadata %>% 
  filter(sample_accession %in% c(droplet_samples, well_samples)) %>% 
  select(study_accession, Platform, sample_accession, Covariate, integration_group) %>% 
  unique() %>% 
  mutate(study_accession = paste0(study_accession, '__', Platform, '__', Covariate, '__', integration_group)) %>%
  mutate(tech = case_when(sample_accession %in% droplet_samples ~ 'droplet',
                          TRUE ~ 'well')) %>% 
  arrange(study_accession)

# drop SRP161678 for now as it's giving weird errors right now
# study_sample <- study_sample %>% filter(!grepl('SRP161678', study_accession))

# merge droplet and well data into one list of sparse matrices by study_sample$study_accession
# option to downsample to no more than 10,000 per study
study_data <- list()
# trycatch for SCTransform
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
for (i in unique(study_sample %>% pull(study_accession))){
  print(i)
  samples <- study_sample %>% filter(study_accession == i) %>% pull(sample_accession) 
  tech <- study_sample %>% filter(study_accession == i) %>% pull(tech) %>% head(1)
  if (tech == 'droplet'){
    study_data[[i]] <- do.call(cbind, sc_data[samples])
  }
  else {
    # well
    # remove samples that aren't in count (which means they failed upstream QC)
    samples <- samples[samples %in% well_samples]
    study_data[[i]] <- count[,samples]
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
  study_data[[i]] <- CreateSeuratObject(study_data[[i]], project = i)
  # calc percentage mito genes
  study_data[[i]][["percent.mt"]] <- PercentageFeatureSet(study_data[[i]], pattern = "^MT-")
  # remove cells with > 10% mito genes
  study_data[[i]] <- subset(study_data[[i]], subset = percent.mt < 10)
  study_data[[i]] <- trySCTransform(study_data[[i]])
  #study_data[[i]] <- RunPCA(study_data[[i]])
}
save(study_data, file = paste0(stamp, '__study_data__emergency_01.Rdata'), compress = FALSE)
# identify sets with less than 200 cells, which will fail integration
low_n <- c()
for (i in names(study_data)){
  if (ncol(study_data[[i]]) < 200){
    low_n <- c(low_n, i)
  }
}

if (length(low_n) > 0){
  study_data[low_n] <- NULL
}

# split again, by "integration group"
# for example, will integrate embryonic and postnatal separately for mouse
groups <- names(study_data) %>% str_split(., '__') %>% map(., 4) %>% unlist() %>% unique()
if (length(groups) == 1) {
  study_data_features <- SelectIntegrationFeatures(object.list = study_data, nfeatures = 3000, verbose = FALSE)
  study_data <- PrepSCTIntegration(object.list = study_data, anchor.features = study_data_features, verbose = FALSE)
  
  # have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
  study_data <- lapply(X = study_data, FUN = function(x) {
    x <- RunPCA(x, features = study_data_features, verbose = FALSE)
  })
  save(study_data_features, study_data, file = paste0(stamp, '__study_data__emergency_02.Rdata'), compress = FALSE)
  
  # try clear some memory
  gc()
  anchors <- FindIntegrationAnchors(object.list = study_data, 
                                    normalization.method = 'SCT', 
                                    #reference = grep('SRP158081__10xv2', names(study_data)),
                                    scale = FALSE, 
                                    anchor.features = study_data_features, 
                                    reduction = "cca")
  
} else {
  anchors <- list()
  for (i in groups){
    study_names <- grep(i, names(study_data), value = T)
    study_data_subset <- study_data[study_names]
    study_data_features <- SelectIntegrationFeatures(object.list = study_data_subset, nfeatures = 3000, verbose = FALSE)
    study_data_subset <- PrepSCTIntegration(object.list = study_data_subset, anchor.features = study_data_features, verbose = FALSE)
    
    # have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
    study_data_subset <- lapply(X = study_data_subset, FUN = function(x) {
      x <- RunPCA(x, features = study_data_features, verbose = FALSE)
    })
    
    
    # try clear some memory
    gc()
    print(paste('Start Anchor Finding with', i ))
    anchors[[i]] <- FindIntegrationAnchors(object.list = study_data_subset, 
                                               normalization.method = 'SCT', 
                                               scale = FALSE, 
                                               anchor.features = study_data_features, 
                                               reduction = "cca")
    
    
  }
}

save(anchors, file = args[1], compress = FALSE)


