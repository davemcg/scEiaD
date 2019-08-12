library(Matrix)
library(liger)
library(tidyverse)
setwd('/data/mcgaugheyd/projects/nei/mcgaughey/massive_integrated_eye_scRNA/')
tx <- read_tsv('references/gencode.vM22.metadata.MGI_tx_mapping.tsv', col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

rdata_files <- c('quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata' ,'quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

metadata <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv')
# load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')
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
}



# liger time
ligerex = createLiger(study_data)
ligerex = normalize(ligerex)
ligerex = selectGenes(ligerex, var.thresh = 0.1, do.plot = FALSE)
ligerex <- scaleNotCenter(ligerex)

# suggestK(ligerex)
ligerex <- optimizeALS(ligerex, k = 25)
ligerex <- quantileAlignSNF(ligerex)

ligerex <- runUMAP(ligerex)
liger_seurat <- ligerToSeurat(ligerex, by.dataset = TRUE)
liger_seurat <- NormalizeData(liger_seurat)
liger_seurat <- RunPCA(liger_seurat, npcs = 100)
liger_seurat <- RunUMAP(liger_seurat, dims=1:20)

save(liger_seurat, ligerex, file = 'liger_ligerex_and_ligerSeurat_2019_08_12.Rdata')
