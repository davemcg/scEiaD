# scEiaD

 Harmonize publicly **available** scRNA-seq datasets for the eye
 
 The code base that underlies the data for plae.nei.nih.gov

# Version 0.92

## Selected parameters for integration

  - scANVI (scVI version 0.13)
  - 4000 HVG
  - 5 epochs
  - 15 latent dimensions

## Yaml files

 - config.yaml for SnakeQUANT and SnakePOP
 - configSCEIAD.yaml for SnakeSCEIAD

 # Numbers
 - 44 datasets
 - 4 species (human, mouse, macaque, chicken)
 - 1,136,041 cells
 - ~60 different cell types
 
 <img width="329" alt="image" src="https://user-images.githubusercontent.com/10225430/189389690-13381192-79fa-48df-945f-8a439bcfcb09.png">

 

# Data
  - [scEiaD database for PLAE](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/MOARTABLES__anthology_limmaFALSE___4000-counts-universe-study_accession-scANVIprojection-15-5-0.1-50-20__pointRelease01.sqlite.gz)
  - [Seurat (v3)](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/scEiaD_all_seurat_v3.Rdata)
  - [anndata (h5ad)](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/scEiaD_all_anndata.h5ad)
  - [Cell Level Metadata](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/metadata_filter.tsv.gz)
  - [Cell Level Metadata in FST format used by web app](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/meta_filter.fst)
  - [Counts (Rdata sparse matrix)](http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/counts.Rdata)

# Moar data?
Visit plae.nei.nih.gov -> Data for differential gene expression tables, study level counts, and more. 

# How do I add my own data to scEiaD?
It is possible! Even better, you don't have to ask me! Or tell me! 

[Interactive colab notebook](https://colab.research.google.com/github/davemcg/scEiaD/blob/master/colab/Query_scEiaD_with_scVI.ipynb)
