# scEiaD

 Harmonize all publicly **available** scRNA-seq datasets for retina
 
 # Numbers
 - 34 datasets
 - 3 species (human, mouse, macaque)
 - 1.2e6 cells go in, 7.7e5 come out
 - 30 different cell types
 
 # Methods
  - Curate data, curate cell type labels
  - Re-process all with [kallisto-bustools](https://www.kallistobus.tools)
  - Grid search and [scPOP](https://github.com/vinay-swamy/scPOP) to ID best integration method (winner: scVI)
  - Grid search again with scVI only to ID best params (5000 HVG, 8 latent dims)
    - Build scVI model on human data, then use scVI/scArches [query mode](https://docs.scvi-tools.org/en/0.9.0/user_guide/notebooks/scarches_scvi_tools.html) to add mouse and macaque data.
  - UMAP
  - xgboost cell type prediction to label unlabelled cells

# Data
  - [Seurat (v3)](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/scEiaD_all_seurat_v3.Rdata)
  - [anndata (h5ad)](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/scEiaD_all_anndata.h5ad)
  - [Cell Level Metadata](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/metadata_filter.tsv.gz)
  - [Counts (Rdata sparse matrix)](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/counts.Rdata)
  
# How do I add my own data to scEiaD?
It is possible! Even better, you don't have to ask me! Or tell me! 

[Interactive colab notebook](https://colab.research.google.com/github/davemcg/scEiaD/blob/master/colab/Query_scEiaD_with_scVI.ipynb)
