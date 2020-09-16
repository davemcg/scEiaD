args = commandArgs(trailingOnly=TRUE)
working_dir = args[1]
git_dir = args[2]

library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(glue)
setwd(working_dir)
load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('x', env)]
  stopifnot(length(var) == 1)
  return(get(var))
}



maca_mf_matrix <- load_rdata('pipeline_data/clean_quant/Macaca_fascicularis/mf-macaca_mulatta_full_sparse_matrix.Rdata')
maca_hs_matrix <- load_rdata('pipeline_data/clean_quant/Macaca_fascicularis/hs-homo_sapiens_full_sparse_matrix.Rdata')




gene_id_converter <- read_tsv('references/ensembl_biomart_human2mouse_macaque.tsv', skip = 1,
                              col_names= c('hs_gene_id','hs_gene_id_v', 'mm_gene_id', 'mf_gene_id',
                                           'hs_gene_name', 'mf_gene_name', 'mm_gene_name')) %>% 
  select(-hs_gene_id_v)

#### This first part is to blend the macaque and human counts 


## first calculate the total counts across each build
maca_mf_rowSums <- tibble(mf_gene_id = rownames(maca_mf_matrix),  mf_total = rowSums(maca_mf_matrix))
maca_hs_rowSums <- tibble(hs_gene_id = rownames(maca_hs_matrix),  hs_total = rowSums(maca_hs_matrix))

## merge counts together
joined <- gene_id_converter %>% select(hs_gene_id, mf_gene_id) %>% distinct %>% 
  left_join(maca_mf_rowSums) %>% left_join(maca_hs_rowSums)
## pick a minimum threshold the total counts must be in order to be considered for blending; this is stop genes that have 
## a few counts from being used; the threshold I've picked is counts >  the first quantile of nonzero gene expression of 
## the macaque annotation
min_hs_exp <- maca_mf_rowSums%>% filter(mf_total > 0) %>% pull(mf_total) %>% quantile(.25)


## next, pick macaque genes that have greater expression than human genes;
mf_genes_id_greater <- joined %>% filter(hs_total <= mf_total) %>% pull(mf_gene_id)
## check to see whether the better gene also has a mapping to a human gene id
mf_genes_in_hs <- gene_id_converter %>% filter(mf_gene_id %in% mf_genes_id_greater) %>% select(hs_gene_id, mf_gene_id) %>% pull(mf_gene_id)
hs_genes_outright_better <- joined %>% filter(hs_total > mf_total, hs_total > min_hs_exp) %>% pull(hs_gene_id)

## for macaque genes that had a higher expression than human genes, but had no mapping, repick human genes again using the expresion threhold
hs_genes_mf_missing <- joined %>% filter(mf_gene_id %in% mf_genes_id_greater, !mf_gene_id %in% mf_genes_in_hs, hs_total > min_hs_exp) %>% pull(hs_gene_id)
hs_genes <- c(hs_genes_outright_better, hs_genes_mf_missing)

## convert macaque gene ids to human gene ids and fix rownames accordingly
hs_to_mf <- gene_id_converter %>% 
  filter(mf_gene_id %in% mf_genes_in_hs) %>% 
  select(hs_gene_id, mf_gene_id) %>% 
  distinct %>% 
  filter(!duplicated(mf_gene_id), mf_gene_id %in% rownames(maca_mf_matrix), 
         )# multiple human gene ids map to the same macaque id, so remove duplicates

## multiple macaque genes also map to the same human gene_id, so sum those together 
maca_mf_matrix_hs_genes <- maca_mf_matrix[hs_to_mf$mf_gene_id, ]
maca_mf_matrix_hs_genes <-  aggregate.Matrix(maca_mf_matrix_hs_genes, groupings = hs_to_mf$hs_gene_id, fun = 'sum')

## keep only the human genes that meet criteria 
maca_hs_matrix_hs_genes <- maca_hs_matrix[rownames(maca_hs_matrix)%in% hs_genes , ]

nrow(maca_hs_matrix_hs_genes) + nrow(maca_mf_matrix_hs_genes)

## since the columns dont match, first identify the set diffs and intersections 
mf_specific_cells <- colnames(maca_mf_matrix_hs_genes)[!colnames(maca_mf_matrix_hs_genes) %in% colnames(maca_hs_matrix_hs_genes)]
hs_specific_cells <- colnames(maca_hs_matrix_hs_genes)[!colnames(maca_hs_matrix_hs_genes) %in% colnames(maca_mf_matrix_hs_genes)]
shared_cells <- colnames(maca_mf_matrix_hs_genes)[colnames(maca_mf_matrix_hs_genes) %in% colnames(maca_hs_matrix_hs_genes)]

## select common cells, then row bind; this lets us have all the genes we need; then RowMerge to fill in missing columns
all_cells_mf_hs_matrix <- rbind(maca_mf_matrix_hs_genes[,shared_cells], maca_hs_matrix_hs_genes[,shared_cells]) %>% 
  RowMergeSparseMatrices(maca_mf_matrix_hs_genes[,mf_specific_cells]) %>% 
  RowMergeSparseMatrices(maca_hs_matrix_hs_genes[,hs_specific_cells])

all_cells_mf_hs_matrix <-  all_cells_mf_hs_matrix[rowSums(all_cells_mf_hs_matrix) >0, ]

all_cells_macaque_hs_ids <- all_cells_mf_hs_matrix
save(all_cells_macaque_hs_ids, file ='pipeline_data/clean_quant/Macaca_fascicularis/all_macaque_full_sparse_matrix.Rdata')

# rebuild cell info.
#### Done with macaque, now lets merge across all species
## free up some memory
gdata::keep(all_cells_mf_hs_matrix, gene_id_converter, joined, load_rdata,  sure = T)
## First, convert mouse ID's to human ids 
homo_hs_matrix <-load_rdata('pipeline_data/clean_quant/Homo_sapiens/hs-homo_sapiens_full_sparse_matrix.Rdata')
mus_mm_matrix <- load_rdata('pipeline_data/clean_quant/Mus_musculus/mm-mus_musculus_full_sparse_matrix.Rdata')

all_shared_gene_ids_hs_mm <- gene_id_converter %>% 
  filter(hs_gene_id %in% rownames(all_cells_mf_hs_matrix), 
         !is.na(mm_gene_id))%>%
  filter(!duplicated(mm_gene_id), # remove mouse genes that map to the same gene ID
         mm_gene_id %in% rownames(mus_mm_matrix)) %>% 
  select(hs_gene_id, mm_gene_id) %>% distinct

## merge everything together 
mus_mm_matrix_cg <- mus_mm_matrix[all_shared_gene_ids_hs_mm$mm_gene_id, ]

mus_mm_matrix_cg <- aggregate.Matrix(mus_mm_matrix_cg, all_shared_gene_ids_hs_mm$hs_gene_id, fun='sum')

homo_hs_matrix_cg <- homo_hs_matrix[rownames(mus_mm_matrix_cg), ]
rm(mus_mm_matrix, homo_hs_matrix)# free up more mem 
all_cells_all_species_matrix <-  RowMergeSparseMatrices(homo_hs_matrix_cg, mus_mm_matrix_cg) %>%  
  RowMergeSparseMatrices(all_cells_mf_hs_matrix)

metadata <- read_tsv(glue('{git_dir}/data/sample_run_layout_organism_tech.tsv'))

all_cell_info <- colnames(all_cells_all_species_matrix) %>% enframe() %>% 
  mutate(sample_accession = str_extract(value, '(ERS|SRS|iPSC_RPE_scRNA_)\\d+')) %>% 
  left_join(metadata %>% select(-run_accession) %>% unique()) %>% 
  data.frame() %>% 
  mutate(batch = paste(study_accession, Platform, Covariate, sep = '_'),
         batch2 = paste(study_accession, Covariate, sep = '_'),
         batch3 = paste(Platform, Covariate, sep = '_')) %>% 
  select(-name )

save(all_cells_all_species_matrix, file = 'pipeline_data/clean_quant/all_species_full_sparse_matrix.Rdata', compress = F)
write_tsv(all_cell_info, path  = 'pipeline_data/cell_info/all_cell_info.tsv')

                                  