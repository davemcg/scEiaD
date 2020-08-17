args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(RETICULATE_PYTHON = "/data/mcgaugheyd/conda/envs/sceasy/bin/python")
library(sceasy)
library(reticulate)
#loompy <- reticulate::import('loompy')
library(dplyr)

# Rsrcript ~/git/massive_integrated_eye_scRNA/src/seurat_to_h5ad_core.R seurat_obj_file seurat_obj_name output_h5ad

load(args[1])
out <- get(args[2])
sceasy::convertFormat(out, from="seurat", to="anndata",
                       outFile=args[3])
