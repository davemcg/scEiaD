rdata_files <- c('quant/SRS866911/genecount/matrix_scTransform.Rdata','quant/SRS866908/genecount/matrix_scTransform.Rdata','quant/SRS1467254/genecount/matrix_scTransform.Rdata','quant/SRS3971245/genecount/matrix_scTransform.Rdata','quant/SRS3971245/genecount/matrix_scTransform.Rdata','quant/SRS3971245/genecount/matrix_scTransform.Rdata','quant/SRS3971245/genecount/matrix_scTransform.Rdata','quant/SRS3971246/genecount/matrix_scTransform.Rdata','quant/SRS3971246/genecount/matrix_scTransform.Rdata','quant/SRS3971246/genecount/matrix_scTransform.Rdata','quant/SRS3971246/genecount/matrix_scTransform.Rdata','quant/SRS4363764/genecount/matrix_scTransform.Rdata','quant/SRS1467251/genecount/matrix_scTransform.Rdata','quant/SRS1467253/genecount/matrix_scTransform.Rdata','quant/SRS3674976/genecount/matrix_scTransform.Rdata','quant/SRS3674982/genecount/matrix_scTransform.Rdata','quant/SRS3674983/genecount/matrix_scTransform.Rdata','quant/SRS4363765/genecount/matrix_scTransform.Rdata','quant/SRS3674974/genecount/matrix_scTransform.Rdata','quant/SRS3674975/genecount/matrix_scTransform.Rdata','quant/SRS3674985/genecount/matrix_scTransform.Rdata','quant/SRS1467249/genecount/matrix_scTransform.Rdata','quant/SRS3674980/genecount/matrix_scTransform.Rdata','quant/SRS3971244/genecount/matrix_scTransform.Rdata','quant/SRS3971244/genecount/matrix_scTransform.Rdata','quant/SRS3971244/genecount/matrix_scTransform.Rdata','quant/SRS3971244/genecount/matrix_scTransform.Rdata','quant/SRS4363763/genecount/matrix_scTransform.Rdata','quant/SRS4386076/genecount/matrix_scTransform.Rdata','quant/SRS1467250/genecount/matrix_scTransform.Rdata','quant/SRS3674978/genecount/matrix_scTransform.Rdata','quant/SRS3674977/genecount/matrix_scTransform.Rdata','quant/SRS3674988/genecount/matrix_scTransform.Rdata','quant/SRS3674979/genecount/matrix_scTransform.Rdata','quant/SRS3674981/genecount/matrix_scTransform.Rdata','quant/SRS3674984/genecount/matrix_scTransform.Rdata','quant/SRS4386075/genecount/matrix_scTransform.Rdata','quant/SRS1467252/genecount/matrix_scTransform.Rdata','quant/SRS3674987/genecount/matrix_scTransform.Rdata','quant/SRS866912/genecount/matrix_scTransform.Rdata','quant/SRS866910/genecount/matrix_scTransform.Rdata','quant/SRS866909/genecount/matrix_scTransform.Rdata','quant/SRS866907/genecount/matrix_scTransform.Rdata','quant/SRS4363762/genecount/matrix_scTransform.Rdata','quant/SRS866906/genecount/matrix_scTransform.Rdata')


library(Matrix)
library(liger)
library(tidyverse)
setwd('/data/mcgaugheyd/projects/nei/mcgaughey/massive_integrated_eye_scRNA/')
tx <- read_tsv('references/gencode.vM22.metadata.MGI_tx_mapping.tsv', col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')

rdata_files <- c('quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata' ,'quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

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
  #row.names(res_matrix) <- gsub('-','_', row.names(res_matrix))
  colnames(res_matrix) <- make.unique(colnames(res_matrix))
  colnames(res_matrix) <- paste0(colnames(res_matrix), "_", sample_accession)
  row.names(res_matrix) <- make.unique(row.names(res_matrix)) %>% toupper()
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
}

# liger time
ligerex = createLiger(study_data)
ligerex = normalize(ligerex)
ligerex = selectGenes(ligerex, var.thresh = 0.1)
ligerex <- scaleNotCenter(ligerex)

# suggestK(ligerex)
ligerex <- optimizeALS(ligerex, k = 25)
ligerex <- quantileAlignSNF(ligerex)

ligerex <- runUMAP(ligerex)
liger_seurat <- ligerToSeurat(ligerex, by.dataset = TRUE)

save(liger_seurat, ligerex, file = 'liger_ligerex_and_ligerSeurat.Rdata')







sc_data <- list()
for (i in seq(1,length(rdata_files))){v
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  colnames(res_matrix) <- make.unique(colnames(res_matrix))
  colnames(res_matrix) <- paste0(colnames(res_matrix), "_", sample_accession)
  row.names(res_matrix) <- make.unique(row.names(res_matrix))
  sc_data[[sample_accession]] <- res_matrix
}


# seurat
sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  sc_data[[sample_accession]] <- CreateSeuratObject(res_matrix, project = sample_accession)
  sc_data[[sample_accession]] <- SCTransform(sc_data[[sample_accession]], variable.features.n = 3000)
  sc_data[[sample_accession]] <- RunPCA(sc_data[[sample_accession]])
}









sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  colnames(res_matrix) <- make.unique(colnames(res_matrix))
  row.names(res_matrix) <- make.unique(row.names(res_matrix))
  sc_data[[sample_accession]] <- res_matrix
}
ligerex = createLiger(sc_data)
sc_data <- list()
for (i in seq(1,length(rdata_files))){
  file = rdata_files[i]
  sample_accession = str_extract(file, '(SRS|iPSC_RPE_scRNA_)\\d+')
  load(file)
  colnames(res_matrix) <- make.unique(colnames(res_matrix))
  colnames(res_matrix) <- paste0(colnames(res_matrix), "_", sample_accession)
  row.names(res_matrix) <- make.unique(row.names(res_matrix))
  sc_data[[sample_accession]] <- res_matrix
}

