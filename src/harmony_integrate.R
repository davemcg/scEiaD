rdata_files <- c('quant/SRS866911/genecount/matrix.Rdata','quant/SRS866908/genecount/matrix.Rdata','quant/SRS1467254/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971245/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS3971246/genecount/matrix.Rdata','quant/SRS4363764/genecount/matrix.Rdata','quant/SRS1467251/genecount/matrix.Rdata','quant/SRS1467253/genecount/matrix.Rdata','quant/SRS3674976/genecount/matrix.Rdata','quant/SRS3674982/genecount/matrix.Rdata','quant/SRS3674983/genecount/matrix.Rdata','quant/SRS4363765/genecount/matrix.Rdata','quant/SRS3674974/genecount/matrix.Rdata','quant/SRS3674975/genecount/matrix.Rdata','quant/SRS3674985/genecount/matrix.Rdata','quant/SRS1467249/genecount/matrix.Rdata','quant/SRS3674980/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS3971244/genecount/matrix.Rdata','quant/SRS4363763/genecount/matrix.Rdata','quant/SRS4386076/genecount/matrix.Rdata','quant/SRS1467250/genecount/matrix.Rdata','quant/SRS3674978/genecount/matrix.Rdata','quant/SRS3674977/genecount/matrix.Rdata','quant/SRS3674988/genecount/matrix.Rdata','quant/SRS3674979/genecount/matrix.Rdata','quant/SRS3674981/genecount/matrix.Rdata','quant/SRS3674984/genecount/matrix.Rdata','quant/SRS4386075/genecount/matrix.Rdata','quant/SRS1467252/genecount/matrix.Rdata','quant/SRS3674987/genecount/matrix.Rdata','quant/SRS866912/genecount/matrix.Rdata','quant/SRS866910/genecount/matrix.Rdata','quant/SRS866909/genecount/matrix.Rdata','quant/SRS866907/genecount/matrix.Rdata','quant/SRS4363762/genecount/matrix.Rdata','quant/SRS866906/genecount/matrix.Rdata')

setwd('/data/mcgaugheyd/projects/nei/mcgaughey/massive_integrated_eye_scRNA')
library(Seurat)
library(harmony)
library(tidyverse)

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

one_giant_matrix <- do.call(cbind, sc_data)
tx <- read_tsv('references/gencode.vM22.metadata.MGI_tx_mapping.tsv', col_names = FALSE) %>% select(2,3) %>% unique()
colnames(tx) <- c('id', 'gene')
row.names(one_giant_matrix) <- row.names(one_giant_matrix) %>% 
  enframe(value = 'id') %>% 
  left_join(., tx, by = 'id') %>% 
  pull(gene)

mouse_retina <- CreateSeuratObject(counts = one_giant_matrix, project = "Mouse Retina", min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = mouse_retina@var.genes, npcs = 20, verbose = FALSE)

mouse_retina@meta.data$SRS <- srs <- str_extract(row.names(mouse_retina@meta.data), "SRS\\d+")

load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')
mouse_retina@meta.data$study_accession <- left_join(mouse_retina@meta.data$SRS %>% enframe(value = 'sample_accession'), sra_metadata_extended %>% select(study_accession, sample_accession) %>% unique()) %>% pull(study_accession)

# https://stackoverflow.com/questions/21382681/kmeans-quick-transfer-stage-steps-exceeded-maximum
gc()
mouse_retina <- mouse_retina %>% RunHarmony("study_accession", max.iter.harmony = 20, epsilon.harmony = -Inf, verbose = TRUE)

mouse_retina <- mouse_retina %>% RunUMAP(reduction = "harmony", dims = 1:20)

save(mouse_retina, file = 'mouse_retina__harmony3.Rdata')