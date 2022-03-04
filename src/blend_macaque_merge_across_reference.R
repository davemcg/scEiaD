#args <- c(getwd(), getwd(), 'references/samplename_patterns.txt')
args = commandArgs(trailingOnly=TRUE)
working_dir = args[1]
git_dir = args[2]
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(glue)
#metadata <- read_tsv(glue('{git_dir}/data/sample_run_layout_organism_tech_biosample_organ_2021_06_07.tsv'))
metadata <- read_tsv(args[3])
setwd(working_dir)
patterns <- scan(args[4], what = character(), sep='\n') %>% paste0(collapse = '|')
genes_to_retain_from_giga_HVG <- read_csv('~/git/scEiaD/data/scEiaD_2021_03_17_GIGA_HVG.csv', col_names = FALSE) %>% pull(1)
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

# ortholog info
ensembl_ortho <- read_tsv(glue("{git_dir}/data/ensembl105_orthologues_hs_mm_mf_gg_zf.tsv.gz"))
colnames(ensembl_ortho) <- c('hs_gene_id', 'hs_gene_name', 'gg_gene_id', 'gg_gene_name','mm_gene_id','mm_gene_name','mf_gene_id','mf_gene_name','Source','zf_gene_id','zf_gene_name')
ensembl_ortho <- ensembl_ortho %>% select(-Source)


## first calculate the total counts across each build
maca_mf_rowSums <- tibble(mf_gene_id = rownames(maca_mf_matrix),  mf_total = rowSums(maca_mf_matrix))
maca_hs_rowSums <- tibble(hs_gene_id = rownames(maca_hs_matrix),  hs_total = rowSums(maca_hs_matrix))

## merge counts together
joined <- ensembl_ortho  %>% select(hs_gene_id, mf_gene_id) %>% distinct %>% 
  left_join(maca_mf_rowSums) %>% left_join(maca_hs_rowSums) %>% 
  mutate(mf_total = replace_na(mf_total, 0), 
         hs_total = replace_na(hs_total, 0))
## pick a minimum threshold the total counts must be in order to be considered for blending; this is stop genes that have 
## a few counts from being used; the threshold I've picked is counts >  the first quantile of nonzero gene expression of 
## the macaque annotation
min_hs_exp <- maca_mf_rowSums%>% filter(mf_total > 0) %>% pull(mf_total) %>% quantile(.25)


## next, pick macaque genes that have greater expression than human genes;
mf_genes_id_greater <- joined %>% filter(hs_total < mf_total) %>% pull(mf_gene_id) %>% unique
## check to see whether the better gene also has a mapping to a human gene id
mf_genes_in_hs <- ensembl_ortho %>% filter(mf_gene_id %in% mf_genes_id_greater) %>%  pull(mf_gene_id) %>% unique
hs_genes_outright_better <- joined %>% filter(hs_total > mf_total, hs_total > min_hs_exp) %>% pull(hs_gene_id)
hs_genes_outright_better_mfid <- joined %>% filter(hs_total > mf_total, hs_total > min_hs_exp) %>% pull(mf_gene_id)
### at this point, any genes that haven't been found have a human expession greater than macaque exprssion, with counts less than 9,
### but we still want to keep any protein coding genes that have been found 


## for macaque genes that had a higher expression than human genes, but had no mapping, repick human genes again using the expresion threhold
hs_genes_mf_missing <- joined %>% filter(mf_gene_id %in% mf_genes_id_greater, !mf_gene_id %in% mf_genes_in_hs, hs_total > min_hs_exp) %>% pull(hs_gene_id)
hs_genes <- c(hs_genes_outright_better, hs_genes_mf_missing)

## convert macaque gene ids to human gene ids and fix rownames accordingly
hs_to_mf <- ensembl_ortho  %>% 
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
## not adding back missing human PC's just yet, I'll let david add them back in with the code I sent
## any missing macaque PC genes have essentially 0 counts

save(all_cells_macaque_hs_ids, file ='pipeline_data/clean_quant/Macaca_fascicularis/full_sparse_matrix.Rdata')


## free up some memory
gdata::keep(all_cells_macaque_hs_ids, joined, load_rdata, git_dir, working_dir, hs_to_mf, 
            hs_genes,merge_macaque_references,patterns, genes_to_retain_from_giga_HVG, sure = T)
#args = commandArgs(trailingOnly=TRUE)
metadata <- read_tsv(args[3])

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
row.names(homo_hs_matrix_cg) <- gsub('\\.\\d+', '', row.names(homo_hs_matrix_cg))
save(homo_hs_matrix_cg, file ='pipeline_data/clean_quant/Homo_sapiens/full_sparse_matrix.Rdata' )
homo_hs__keep_genes = rownames(homo_hs_matrix_cg)
rm(homo_hs_matrix_cg)

### now make the merged quantfile
all_shared_gene_ids <- ensembl_ortho  %>% 
  filter(#hs_gene_id %in% rownames(all_cells_macaque_hs_ids),
         hs_gene_id %in% rownames(homo_hs_matrix),
         mm_gene_id %in% rownames(mus_mm_matrix),
         !is.na(mm_gene_id))#%>%
#  filter(!duplicated(mm_gene_id), # remove mouse genes that map to the same gene ID
#         mm_gene_id %in% rownames(mus_mm_matrix)) #%>% 
  

all_shared_gene_ids_hs_mm <- all_shared_gene_ids %>% select(hs_gene_id, mm_gene_id) %>% distinct

# two "alt" situations remain
## one human gene maps to multiple mouse genes 
## one mouse gene maps to multiple human genes

## let us remove those from the match
## missing genes in a species can be added later for diff testing / viz / etc

## mouse genes that multi-match human genes
mm_multimap <- all_shared_gene_ids_hs_mm %>% group_by(mm_gene_id) %>% summarise(Count = n()) %>% filter(Count > 1)
## human
hs_multimap <- all_shared_gene_ids_hs_mm %>% group_by(hs_gene_id) %>% summarise(Count = n()) %>% filter(Count > 1)

all_shared_gene_ids_hs_mm_one_to_one <- all_shared_gene_ids_hs_mm %>% filter(!mm_gene_id %in% mm_multimap$mm_gene_id) %>%
																	  filter(!hs_gene_id %in% hs_multimap$hs_gene_id)
## IF (big IF) there is a 1 to 1 match between human  name and mouse name then use that 
one_to_ones <- all_shared_gene_ids %>% 
				filter(
					(mm_gene_id %in% 
					mm_multimap$mm_gene_id) |
					(hs_gene_id %in%
				    hs_multimap$hs_gene_id)) %>% 
				filter(hs_gene_name == toupper(mm_gene_name)) %>% 
				select(hs_gene_id, mm_gene_id) %>% 
				distinct()

all_shared_gene_ids_hs_mm_one_to_one <- bind_rows(all_shared_gene_ids_hs_mm_one_to_one, one_to_ones)

# again check for oen to many
again_check <- bind_rows(one_to_ones %>% group_by(mm_gene_id) %>% summarise(Count = n()) %>% filter(Count > 1), one_to_ones %>% group_by(hs_gene_id) %>% summarise(Count = n()) %>% filter(Count > 1))

all_shared_gene_ids_hs_mm_one_to_one <- all_shared_gene_ids_hs_mm_one_to_one %>% filter(!hs_gene_id %in% again_check$hs_gene_id) %>% filter(!mm_gene_id %in% again_check$mm_gene_id)

## cut down to shared one to one genes
mus_mm_matrix_cg <- mus_mm_matrix[all_shared_gene_ids_hs_mm_one_to_one$mm_gene_id, ]
row.names(mus_mm_matrix_cg) <- all_shared_gene_ids_hs_mm_one_to_one$hs_gene_id
# fix row names for homo gene names (remove .\\d+ endings)
row.names(homo_hs_matrix) <- str_replace_all(row.names(homo_hs_matrix),'\\.\\d+', '')
homo_hs_matrix_cg <- homo_hs_matrix[rownames(mus_mm_matrix_cg), ]

# cut down macaque matrix to just gene IDs in mus_mm_matrix_cg (which is now using ENSG huma nids)
maca_all_matrix_cg = all_cells_macaque_hs_ids[rownames(mus_mm_matrix_cg)[rownames(mus_mm_matrix_cg) %in% row.names(all_cells_macaque_hs_ids)],  ]

rm(mus_mm_matrix, homo_hs_matrix)# free up more mem 

## add extra species
## right now just chick (gallus gallus)
chick_gg_matrix <-load_rdata('pipeline_data/clean_quant/Gallus_gallus/gg-gallus_gallus_full_sparse_matrix.Rdata')
save(chick_gg_matrix, file ='pipeline_data/clean_quant/Gallus_gallus/full_sparse_matrix.Rdata')
print('extra species time')
extra_species_metadata <- tibble(species = 'Gallus gallus', 
                                 file = 'pipeline_data/clean_quant/Gallus_gallus/gg-gallus_gallus_full_sparse_matrix.Rdata', 
                                 prefix = 'gg_')

process_extra_species <- function(species, file, prefix, converter_table){
   print('load data')
   species_counts <- load_rdata(file) 
   spec_col <- paste0(prefix, 'gene_id')
   converter_cols <- c('hs_gene_id', spec_col )
   print('converter table')
   converter_table <- converter_table %>%
	 select(hs_gene_id, all_of(spec_col)) %>% unique() %>% 
     filter(!is.na(.[,spec_col]),
            .[[spec_col]] %in% rownames(species_counts)
            ) 
   # remove one to many from spec to hs
   one_to_many <- converter_table %>% group_by_at(vars(all_of(spec_col))) %>% summarise(Count = n()) %>% filter(Count > 1)
   converter_table %>% filter(!.data[[spec_col]] %in% (one_to_many %>% pull(1)) )
	print('filter')
   species_counts_filtered <- species_counts[converter_table[[spec_col]],]
   rownames(species_counts_filtered) <- converter_table$hs_gene_id
   # in cases where multiple chick genes match one human gene, sum the chick gene counts
   # (aggregated by human gene id)
   species_counts_filtered <- aggregate.Matrix(species_counts_filtered, groupings = rownames(species_counts_filtered), fun = 'sum')
   return(species_counts_filtered)
}
extra_species <- lapply(split(extra_species_metadata, 1:nrow(extra_species_metadata)), 
       function(x) process_extra_species(x$species, 
                                         x$file, 
                                         x$prefix, 
                                         all_shared_gene_ids)
       )


chick_all_matrix_cg <- extra_species[[1]] 
save(homo_hs_matrix_cg, mus_mm_matrix_cg, maca_all_matrix_cg, chick_all_matrix_cg, file = 'pipeline_data/clean_quant/gene_name_norm.mat.Rdata')

# merge together (THIS IS WHERE YOU GET THE CHOLMOD / vec ISSUES)
all_species_matrices <- c(list(homo_hs_matrix_cg,mus_mm_matrix_cg, maca_all_matrix_cg), 
                          chick_all_matrix_cg)

all_cells_all_species_matrix <-  reduce(all_species_matrices, RowMergeSparseMatrices)


all_cell_info <- colnames(all_cells_all_species_matrix) %>% enframe() %>% 
  mutate(sample_accession = str_extract(value, glue('({patterns})\\d+') )) %>% 
  left_join(metadata %>% select(-run_accession) %>% unique()) %>% 
  data.frame() %>% 
  mutate(batch = paste(study_accession, Platform, Covariate, sep = '_'),
         batch2 = paste(study_accession, Covariate, sep = '_'),
         batch3 = paste(Platform, Covariate, sep = '_')) %>% 
  select(-name )

# grab gene names
all_species_full_sparse_genes <- row.names(all_cells_all_species_matrix)
save(all_species_full_sparse_genes, all_species_full_sparse_remainder_genes, file = 'pipeline_data/clean_quant/gene_names.Rdata')
print('save big matrix')
save(all_cells_all_species_matrix, file = 'pipeline_data/clean_quant/all_species_full_sparse_matrix.Rdata', compress = F)

print('save all cell info')
write_tsv(all_cell_info, path  = 'pipeline_data/cell_info/all_cell_info.tsv')
gene_id_converter %>% select(hs_gene_id, hs_gene_name) %>% distinct %>% write_tsv('references/ENSG2gene_name.tsv.gz')

stats_files <- list.files('pipeline_data/clean_quant', pattern = 'stats.tsv', recursive= T, full.names=T) %>% 
  .[!grepl('droplet_quant_stats.tsv', .)]
studies = str_split(stats_files, '/') %>% sapply(function(x)x[3])
all_stats <- lapply(seq_along(stats_files), function(i) read_tsv(stats_files[i]) %>%  
                      mutate(study_accession = studies[i]) %>% 
                      select(study_accession, everything())) %>% bind_rows 
 

write_tsv(all_stats, 'pipeline_data/clean_quant/droplet_quant_stats.tsv')

