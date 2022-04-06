args <- commandArgs(trailingOnly = TRUE)
library(glue)
conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')

Sys.setenv(RETICULATE_PYTHON = glue("{conda_dir}envs/sceasy/bin/python"))
library(sceasy)
library(reticulate)
#loompy <- reticulate::import('loompy')
library(dplyr)
# Rsrcript ~/git/massive_integrated_eye_scRNA/src/seurat_to_h5ad_core.R seurat_obj_file seurat_obj_name output_h5ad

set.seed(326490)
load(args[1])
out <- get(args[2])

load('pipeline_data/cell_info/cell_info_labelled.Rdata')

nmeta <- out@meta.data %>% as_tibble(rownames = 'value') %>% select(-batch, -Age, -study_accession) %>% left_join(cell_info_labels, by = 'value')
nmeta_df <- data.frame(nmeta)
row.names(nmeta_df) <- nmeta$value

out@meta.data <- nmeta_df

 
if (args[4] == 'make_mini_split_data'){
	ref_samples <- scan(glue('{git_dir}data/human_ref_samples.txt'), what = 'character')
	ref_bc_mini <- scEiaD@meta.data %>% as_tibble(rownames = 'Barcode2') %>% filter(sample_accession %in% ref_samples) %>% group_by(batch) %>% sample_n(2000, replace = TRUE) %>% unique() %>% pull(Barcode2)
	query_bc_mini <- scEiaD@meta.data %>% as_tibble(rownames = 'Barcode2') %>% filter(!sample_accession %in% ref_samples, study_accession != 'SRP131661') %>% group_by(batch) %>% sample_n(500, replace = TRUE) %>% unique() %>% pull(Barcode2)
	out_ref = out[, ref_bc_mini]
	out_query = out[, query_bc_mini]
	sceasy::convertFormat(out_ref, from="seurat", to="anndata",
                       outFile=gsub('.h5ad','_mini_ref.h5ad', args[3]))
	sceasy::convertFormat(out_query, from="seurat", to="anndata",
                       outFile=gsub('.h5ad','_mini_query.h5ad', args[3]))
} 
sceasy::convertFormat(out, from="seurat", to="anndata",
                       outFile=args[3])
