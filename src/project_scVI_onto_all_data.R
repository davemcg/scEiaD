args <- commandArgs(trailingOnly = TRUE)
conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
library(glue)
Sys.setenv(RETICULATE_PYTHON = glue("{conda_dir}/envs/sceasy/bin/python"))
library(sceasy)
library(reticulate)
#loompy <- reticulate::import('loompy')
library(dplyr)
library(Seurat)

# Rsrcript ~/git/massive_integrated_eye_scRNA/src/seurat_to_h5ad_core.R seurat_obj_file seurat_obj_name output_h5ad

# convert seurat object to h5ad for scVI
load(args[1])
existing_seurat_obj <- get(args[5])
model_dir <- args[2]
seurat_obj_out_name <- args[4]
# if args[2] is not None, then new dat is being provided that needs to be projected with 
# the scVI model
# the "existing" data is the data used to train the scVI model
if (args[3] != 'None'){
	load(args[3])
	new_seurat_obj <-  get(args[6])
	h5ad_file_existing = paste0(args[1], 'existing.h5ad')
	
	h5ad_file_new = paste0(args[3], 'new.h5ad')
	scVI_meta_out <- glue("{h5ad_file_new}_meta.csv")
	scVI_ld_out <- glue("{h5ad_file_new}_ld.csv")
	# convert to h5ad for scVI loading
	sceasy::convertFormat(existing_seurat_obj, from="seurat", to="anndata",
                       outFile=h5ad_file_existing)
	sceasy::convertFormat(new_seurat_obj, from="seurat", to="anndata",
                       outFile=h5ad_file_new)

	# now run scVI to load the existing model and output the latent dims
	# for the data loaded here
	scvi_job <- glue("{conda_dir}/envs/scvitools/bin/python {git_dir}/src/project_scVI_onto_all_data.py {h5ad_file_existing} {h5ad_file_new} {scVI_meta_out} {scVI_ld_out} {model_dir}")
	system(scvi_job)
	# load in scVI ld and meta to put scVI in seurat_obj
	meta <- readr::read_csv(scVI_meta_out)
	ld <- readr::read_csv(scVI_ld_out)

	metald <- cbind(meta, ld[,2:ncol(ld)])

    cat_seurat <- MergeSeurat(object1 = existing_seurat_obj, object2 = new_seuerat_obj, 
								add.cell.id1 = "reference", add.cell.id2 = "query")
	ld_ordered <- cat_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% left_join(metald %>% mutate(Barcode =index), by = c('Barcode'))
	scVI_cols <- grep('^\\d+$', colnames(ld_ordered))
	ld_ordered <-  ld_ordered[,scVI_cols] %>% as.matrix()
	row.names(ld_ordered) <- row.names(cat_obj@meta.data)
	colnames(ld_ordered) <- paste0('scVI_',seq(1,ncol(ld_ordered)))
	cat_seurat[["scVI"]] <- CreateDimReducObject(embeddings = ld_ordered, key = "scVI_", assay = DefaultAssay(cat_seurat))

	integrated_obj <- cat_obj
	save(integrated_obj, file = seurat_obj_out_name)
} else {
	integrated_obj <- existing_seurat_obj
	save(integrated_obj, file = seurat_obj_out_name)
}
