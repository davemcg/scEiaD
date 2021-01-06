args = commandArgs(trailingOnly=TRUE)
working_dir = args[1]
git_dir = args[2]
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(glue)
setwd(working_dir)
patterns <- scan(args[3], what = character(), sep='\n') %>% paste0(collapse = '|')
load_rdata <- function(x){
  load(x)
  env <- ls.str()
  var <- env[!grepl('^x$', env)]
  stopifnot(length(var) == 1)
  all_data= get(var)# the mistake was here, I orignally fixing rownames in here, but deleted it by accident
  return(all_data)
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
mf_genes_id_greater <- joined %>% filter(hs_total < mf_total) %>% pull(mf_gene_id) %>% unique
## check to see whether the better gene also has a mapping to a human gene id
mf_genes_in_hs <- gene_id_converter %>% filter(mf_gene_id %in% mf_genes_id_greater) %>%  pull(mf_gene_id) %>% unique
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

merge_macaque_references  = function(maca_mf_matrix,maca_hs_matrix, hs_to_mf, hs_genes){
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
    return(all_cells_mf_hs_matrix)
}

all_cells_macaque_hs_ids <- merge_macaque_references(maca_mf_matrix,maca_hs_matrix, hs_to_mf, hs_genes)
all_cells_macaque_hs_ids <-  all_cells_macaque_hs_ids[rowSums(all_cells_macaque_hs_ids) >0, ]

save(all_cells_macaque_hs_ids, file ='pipeline_data/clean_quant/Macaca_fascicularis/full_sparse_matrix.Rdata')


## free up some memory

gdata::keep(all_cells_macaque_hs_ids, gene_id_converter, joined, load_rdata, git_dir, working_dir, hs_to_mf, 
            hs_genes,merge_macaque_references, sure = T)
intron_maca_mf_matrix <- load_rdata('pipeline_data/clean_quant/Macaca_fascicularis/mf-macaca_mulatta_full_sparse_unspliced_matrix.Rdata')
intron_maca_hs_matrix <- load_rdata('pipeline_data/clean_quant/Macaca_fascicularis/hs-homo_sapiens_full_sparse_unspliced_matrix.Rdata')
all_intron_macaque_data =  merge_macaque_references(intron_maca_mf_matrix,intron_maca_hs_matrix, hs_to_mf, hs_genes)
all_intron_macaque_data = all_intron_macaque_data[rownames(all_cells_macaque_hs_ids), ]
rm(intron_maca_mf_matrix,intron_maca_hs_matrix )
save(all_intron_macaque_data, file ='pipeline_data/clean_quant/Macaca_fascicularis/full_sparse_unspliced_matrix.Rdata')
#### Done with macaque, now lets merge across all species
## First, convert mouse ID's to human ids 
homo_hs_matrix <-load_rdata('pipeline_data/clean_quant/Homo_sapiens/hs-homo_sapiens_full_sparse_matrix.Rdata')
mus_mm_matrix <- load_rdata('pipeline_data/clean_quant/Mus_musculus/mm-mus_musculus_full_sparse_matrix.Rdata')

## save a mouse only quant file 
mus_mm_matrix_hg = mus_mm_matrix[rowSums(mus_mm_matrix) >0, ]
save(mus_mm_matrix_hg, file ='pipeline_data/clean_quant/Mus_musculus/full_sparse_matrix.Rdata' )
mus_mm__keep_genes = rownames(mus_mm_matrix_hg)
rm(mus_mm_matrix_hg)
## save a human only quantfile
homo_hs_matrix_cg = homo_hs_matrix[rowSums(homo_hs_matrix) > 0, ]
save(homo_hs_matrix_cg, file ='pipeline_data/clean_quant/Homo_sapiens/full_sparse_matrix.Rdata' )
homo_hs__keep_genes = rownames(homo_hs_matrix_cg)
rm(homo_hs_matrix_cg)

### now make the merged quantfile
all_shared_gene_ids_hs_mm <- gene_id_converter %>% 
  filter(hs_gene_id %in% rownames(all_cells_macaque_hs_ids), 
         !is.na(mm_gene_id))%>%
  filter(!duplicated(mm_gene_id), # remove mouse genes that map to the same gene ID
         mm_gene_id %in% rownames(mus_mm_matrix)) %>% 
  select(hs_gene_id, mm_gene_id) %>% distinct

## merge everything together 
mus_mm_matrix_cg <- mus_mm_matrix[all_shared_gene_ids_hs_mm$mm_gene_id, ]

mus_mm_matrix_cg <- aggregate.Matrix(mus_mm_matrix_cg, all_shared_gene_ids_hs_mm$hs_gene_id, fun='sum')

homo_hs_matrix_cg <- homo_hs_matrix[rownames(mus_mm_matrix_cg), ]
maca_all_matrix_cg = all_cells_macaque_hs_ids[rownames(mus_mm_matrix_cg),  ]# BUGFIX: - was not removing nonshared genes from macaque
rm(mus_mm_matrix, homo_hs_matrix)# free up more mem 
all_cells_all_species_matrix <-  RowMergeSparseMatrices(homo_hs_matrix_cg, mus_mm_matrix_cg) %>%  
  RowMergeSparseMatrices(maca_all_matrix_cg)

metadata <- read_tsv(glue('{git_dir}/data/sample_run_layout_organism_tech.tsv'))

all_cell_info <- colnames(all_cells_all_species_matrix) %>% enframe() %>% 
  mutate(sample_accession = str_extract(value, glue('({patterns})\\d+') )) %>% 
  left_join(metadata %>% select(-run_accession) %>% unique()) %>% 
  data.frame() %>% 
  mutate(batch = paste(study_accession, Platform, Covariate, sep = '_'),
         batch2 = paste(study_accession, Covariate, sep = '_'),
         batch3 = paste(Platform, Covariate, sep = '_')) %>% 
  select(-name )

save(all_cells_all_species_matrix, file = 'pipeline_data/clean_quant/all_species_full_sparse_matrix.Rdata', compress = F)


write_tsv(all_cell_info, path  = 'pipeline_data/cell_info/all_cell_info.tsv')
gene_id_converter %>% select(hs_gene_id, hs_gene_name) %>% distinct %>% write_tsv('references/ENSG2gene_name.tsv.gz')

## make intron quant for mouse 
rm(all_cells_all_species_matrix)

intron_mus_mm_matrix <- load_rdata('pipeline_data/clean_quant/Mus_musculus/mm-mus_musculus_full_sparse_unspliced_matrix.Rdata')
intron_mus_mm_matrix  = intron_mus_mm_matrix[rownames(intron_mus_mm_matrix)%in%mus_mm__keep_genes, ]
save(intron_mus_mm_matrix, file ='pipeline_data/clean_quant/Mus_musculus/full_sparse_unspliced_matrix.Rdata')

## human intron quant 
intron_homo_hs_matrix <-load_rdata('pipeline_data/clean_quant/Homo_sapiens/hs-homo_sapiens_full_sparse_unspliced_matrix.Rdata')

save(intron_homo_hs_matrix, file ='pipeline_data/clean_quant/Homo_sapiens/full_sparse_unspliced_matrix.Rdata')

### merge all intron quant 
mm_intron_keep = rownames(intron_mus_mm_matrix) %in% all_shared_gene_ids_hs_mm$mm_gene_id
all_shared_gene_ids_hs_mm_intron  =  {tibble(mm_gene_id = rownames(intron_mus_mm_matrix))} %>% 
  inner_join(all_shared_gene_ids_hs_mm)
intron_mus_mm_matrix <- aggregate.Matrix(intron_mus_mm_matrix[mm_intron_keep, ],  
                                         all_shared_gene_ids_hs_mm_intron$hs_gene_id, fun='sum')
intron_homo_hs_matrix <- intron_homo_hs_matrix[rownames(intron_homo_hs_matrix)%in% rownames(intron_mus_mm_matrix), ]
all_intron_macaque_data = all_intron_macaque_data[rownames(intron_mus_mm_matrix), ]
all_intron_data = RowMergeSparseMatrices(intron_homo_hs_matrix, intron_mus_mm_matrix) %>%  
  RowMergeSparseMatrices(all_intron_macaque_data)
save(all_intron_data, file ='pipeline_data/clean_quant/all_species_full_sparse_unspliced_matrix.Rdata')

stats_files <- list.files('pipeline_data/clean_quant', pattern = 'stats.tsv', recursive= T, full.names=T) %>% 
  .[!grepl('droplet_quant_stats.tsv', .)]
studies = str_split(stats_files, '/') %>% sapply(function(x)x[3])
all_stats <- lapply(seq_along(stats_files), function(i) read_tsv(stats_files[i]) %>%  
                      mutate(study_accession = studies[i]) %>% 
                      select(study_accession, everything())) %>% bind_rows 
  

write_tsv(all_stats, 'pipeline_data/clean_quant/droplet_quant_stats.tsv')


