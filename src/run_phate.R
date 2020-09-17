conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
Sys.setenv(RETICULATE_PYTHON = paste0(conda_dir, "/envs/phate/bin/python") )
library(phateR)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])
#load('seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umap.Rdata')

phate_2D <- phate(integrated_obj@reductions$scVI@cell.embeddings, knn = 5, n.jobs=24)
save(phate_2D, file = args[2])



