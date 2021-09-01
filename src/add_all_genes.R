# the merge methods to make the pan-species matrices only keeps genes that can be cross-matched
# between mouse/human/macaque
#
# in no surprise, many genes are dropped
# not a big deal for the integration
# but a problem for the diff testing
# and general db use
# 
# we already know that OPN1MW is lost
# POU4F3 is lost
# so here we grab all (protein-coding) genes
# we are using a pan-tx ref, so i donn't want the expense of keeping every lincRNA/whatevers
conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
args <- commandArgs(trailingOnly = TRUE)
library(glue)
library(tidyverse)
#setwd(path/above/quant)
load('pipeline_data/clean_quant/all_species_full_sparse_matrix.Rdata')
load('pipeline_data/clean_quant/Homo_sapiens/full_sparse_matrix.Rdata')
load('pipeline_data/clean_quant/Mus_musculus/full_sparse_matrix.Rdata')
load('pipeline_data/clean_quant/Macaca_fascicularis/full_sparse_matrix.Rdata')
load('pipeline_data/clean_quant/Gallus_gallus/gg-gallus_gallus_full_sparse_matrix.Rdata')
gg_matrix <- all_data
hs_gtf <- rtracklayer::readGFF('references/gtf/hs-homo_sapiens_anno.gtf.gz') %>% 
  mutate(gene_id = str_remove_all(gene_id, '\\.\\d+$'))
mm_gtf <- rtracklayer::readGFF('references/gtf/mm-mus_musculus_anno.gtf.gz') %>% 
  mutate(gene_id = str_remove_all(gene_id, '\\.\\d+$'))
mf_gtf <- rtracklayer::readGFF('references/gtf/mf-macaca_mulatta_anno.gtf.gz')
gg_gtf <- rtracklayer::readGFF('references/gtf/gg-gallus_gallus_anno.gtf.gz')

conv_tab <- read_tsv(glue('{git_dir}/data/ensembl_biomart_human2mouse_macaque_chick_ZF.tsv.gz'), skip = 1,
                              col_names= c('human_gene_id','verisioned_hid',
                                           'chick_gene_id', 'chick_gene_name',
                                           'macaque_gene_id', 'maca_gene_name',
                                           'mouse_gene_id', 'mouse_gene_name',
                                           'zf_gene_id', 'zf_gene_name',
                                           'gene_name')) 

#conv_tab <-  read_tsv('references/ensembl_biomart_human2mouse_macaque.tsv',
#                      skip = 1,
#                      col_names = c('human_gene_id', 'verisioned_hid','mouse_gene_id', 'macaque_gene_id', 'gene_name', 'mouse_gene_name', 'maca_gene_name')) 
any_common_ctab <- full_join(   
  conv_tab %>% filter(!duplicated(human_gene_id), !duplicated(mouse_gene_id), !is.na(mouse_gene_id)) %>%
    select(human_gene_id, mouse_gene_id),
  conv_tab %>% filter(!duplicated(human_gene_id), !duplicated(macaque_gene_id), !is.na(macaque_gene_id)) %>%
    select(human_gene_id, macaque_gene_id)
) %>%
  full_join(.,
  	conv_tab %>% filter(!duplicated(human_gene_id), !duplicated(chick_gene_id), !is.na(chick_gene_id)) %>%
    select(human_gene_id, chick_gene_id)
  ) %>% inner_join(conv_tab %>% select(gene_name, human_gene_id) %>% distinct,.)

metadata <- read_tsv(glue('{git_dir}data/sample_run_layout_organism_tech.tsv'))
hs_msg_gene_ids <- hs_gtf %>% 
  filter(type == 'gene', !gene_id %in% any_common_ctab$human_gene_id, gene_type == 'protein_coding') %>% 
  select(gene_name, human_gene_id=gene_id) %>% 
  mutate(mouse_gene_id = NA, macaque_gene_id= NA, chick_gene_id = NA)
mm_msg_gene_ids <-  mm_gtf %>% 
  filter(type == 'gene', !gene_id %in% any_common_ctab$mouse_gene_id, gene_type == 'protein_coding') %>%
  select(gene_name, mouse_gene_id=gene_id) %>% 
  mutate(human_gene_id = NA, macaque_gene_id= NA, chick_gene_id = NA, gene_name=toupper(gene_name))
mf_msg_gene_ids <-  mf_gtf %>%
  filter(type == 'gene', !gene_id %in% any_common_ctab$macaque_gene_id, gene_biotype == 'protein_coding') %>%
  select(gene_name, macaque_gene_id=gene_id) %>%
  mutate(human_gene_id = NA, mouse_gene_id= NA, chick_gene_id = NA, gene_name=toupper(gene_name))
gg_msg_gene_ids <- gg_gtf %>%
  filter(type == 'gene', !gene_id %in% any_common_ctab$chick_gene_id, gene_biotype == 'protein_coding') %>%
  select(gene_name, chick_gene_id=gene_id) %>%
  mutate(human_gene_id = NA, mouse_gene_id= NA, macaque_gene_id = NA, gene_name=toupper(gene_name))


all_ids <- bind_rows(any_common_ctab %>% mutate(common=T),
                     hs_msg_gene_ids %>% mutate(common=F),
                     mm_msg_gene_ids %>% mutate(common=F),
					 gg_msg_gene_ids %>% mutate(common=F))

full_map <- all_ids%>% group_by(gene_name) %>% 
  summarise(human_gene_id = first(human_gene_id[!is.na(human_gene_id)]),
            mouse_gene_id = first(mouse_gene_id[!is.na(mouse_gene_id)]),
			maca_gene_id = first(macaque_gene_id[!is.na(macaque_gene_id)]),
			chick_gene_id = first(chick_gene_id[!is.na(chick_gene_id)]),
            common = any(common))
# full_map filtered down to genes in the all_cells_all_species_matrix
merged_map <- full_map %>% filter(human_gene_id %in% rownames(all_cells_all_species_matrix))

# in orig matrix, not in merged matrix
missing_human_quant = homo_hs_matrix_cg[!rownames(homo_hs_matrix_cg) %in% rownames(all_cells_all_species_matrix),]
# is protein coding (present in full_map)
missing_human_quant = missing_human_quant[rownames(missing_human_quant) %in% full_map$human_gene_id,]

# in orig matrix, not in merged matrix
missing_mouse_quant = mus_mm_matrix_hg[!rownames(mus_mm_matrix_hg) %in% merged_map$mouse_gene_id,]
# is protein coding (present in full_map)
missing_mouse_quant = missing_mouse_quant[rownames(missing_mouse_quant) %in% full_map$mouse_gene_id,]
# if gene shared between human anad mouse, then replace ENSMMUS with ENSG (human)
mouse_gene_ids <- rownames(missing_mouse_quant)
new_mouse_gene_ids <- mouse_gene_ids %>% 
		enframe(value = 'mouse_gene_id') %>% 
		select(-name) %>% 
		left_join(full_map) %>%
		mutate(new_id = case_when(is.na(human_gene_id) ~ mouse_gene_id,
									TRUE ~ human_gene_id)) %>% 
		pull(new_id)
rownames(missing_mouse_quant) <- new_mouse_gene_ids

# in orig matrix, not in merged matrix
missing_chick_quant = gg_matrix[!rownames(gg_matrix) %in% merged_map$chick_gene_id,]
# is protein coding (present in full_map)
missing_chick_quant = missing_chick_quant[rownames(missing_chick_quant) %in% full_map$chick_gene_id,]
# if gene shared between human anad mouse, then replace ENSMMUS with ENSG (human)
chick_gene_ids <- rownames(missing_chick_quant)
new_chick_gene_ids <- chick_gene_ids %>%
        enframe(value = 'chick_gene_id') %>%
        select(-name) %>%
        left_join(full_map) %>%
        mutate(new_id = case_when(is.na(human_gene_id) ~ chick_gene_id,
                                    TRUE ~ human_gene_id)) %>%
        pull(new_id)
rownames(missing_chick_quant) <- new_chick_gene_ids



# make macaque missing matrix with the genes in missing_human_quant
# the "all_cells_macaque_hs_ids" has already been labelled with the human gene ids
# it is weird as a hybrid pseudoquant was done where the macaque data was
# quantified against both human and macaque tx-omes
# and human tx-ome was used if the corresponding gene was missing in the macaque
missing_macaque_quant = all_cells_macaque_hs_ids[rownames(missing_human_quant)[ (rownames(missing_human_quant) %in% rownames(all_cells_macaque_hs_ids))], ] 


missing_quant <- Seurat::RowMergeSparseMatrices(missing_human_quant, missing_mouse_quant) %>%
  Seurat::RowMergeSparseMatrices(missing_macaque_quant) %>%
  Seurat::RowMergeSparseMatrices(missing_chick_quant)

load(args[1]) #integrated_obj

counts_integrated_obj <- integrated_obj@assays$RNA@counts
missing_quant_aligned <- missing_quant[, colnames(counts_integrated_obj)]

save(missing_quant, missing_quant_aligned, file = args[2]) #'pipeline_data/clean_quant/missing_quant.Rdata')
save(full_map, file = args[3]) #'pipeline_data/clean_quant/missing_quant_full_map.Rdata')

integrated_obj@assays$RNA@counts <- Matrix::rbind2(counts_integrated_obj, missing_quant_aligned)
# recreate RNA assay to get info right
x <- Seurat::CreateSeuratObject(integrated_obj@assays$RNA@counts, assay = 'RNA')
integrated_obj@assays$RNA <- x@assays$RNA

save(integrated_obj, file = args[4])
