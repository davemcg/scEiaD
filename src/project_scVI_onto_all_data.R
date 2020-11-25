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
seurat_obj <- get(args[4])
h5ad_file = paste0(args[1], '.h5ad')
scVI_meta_out <- glue("{h5ad_file}_meta.csv")
scVI_ld_out <- glue("{h5ad_file}_ld.csv")
model_dir <- args[2]
seurat_obj_out_name <- args[3]
sceasy::convertFormat(seurat_obj, from="seurat", to="anndata",
                       outFile=h5ad_file)

# now run scVI to load the existing model and output the latent dims
# for the data loaded here
scvi_job <- glue("{conda_dir}/envs/scvitools/bin/python {git_dir}/src/project_scVI_onto_all_data.py {h5ad_file} {scVI_meta_out} {scVI_ld_out} {model_dir}")
system(scvi_job)
# load in scVI ld and meta to put scVI in seurat_obj
meta <- readr::read_csv(scVI_meta_out)
ld <- readr::read_csv(scVI_ld_out)

metald <- cbind(meta, ld[,2:ncol(ld)])
ld_ordered <- seurat_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% left_join(metald %>% mutate(Barcode =index), by = c('Barcode'))
scVI_cols <- grep('^\\d+$', colnames(ld_ordered))
ld_ordered <-  ld_ordered[,scVI_cols] %>% as.matrix()
row.names(ld_ordered) <- row.names(seurat_obj@meta.data)
colnames(ld_ordered) <- paste0('scVI_',seq(1,ncol(ld_ordered)))
seurat_obj[["scVI"]] <- CreateDimReducObject(embeddings = ld_ordered, key = "scVI_", assay = DefaultAssay(seurat_obj))

integrated_obj <- seurat_obj
save(integrated_obj, file = seurat_obj_out_name)
