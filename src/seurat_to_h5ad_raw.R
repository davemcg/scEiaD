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
sceasy::convertFormat(out, from="seurat", to="anndata",
                       outFile=args[3])
