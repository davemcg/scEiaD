# massive_integrated_eye_scRNA

Goal: Harmonize all publicly **available** scRNA-seq datasets for eye (mostly retina as of February 2020).

Milestones: 

 - Presented poster Genome Informatics 2019 (November). Managed to get a prototype integration with [fastMNN](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html). 
  - Lior Pachter suggested using [scVI](https://scvi.readthedocs.io/en/stable/) to do the dimension reduction step
  - He also suggested using his/Sina's [kb python](https://github.com/pachterlab/kb_python) to run velocity estimates
- As of Feb 2020 I've swapped over to scVI for the dimensionality reduction and tested a **wide** range of parameters and more rigorously assessed integration performance with ARI, LISI, and some home-spun methods as it was too onerous assessing the UMAP projections by eye. 
- I'm presenting this work at ARVO 2020 (May)! To do stuffs:
  - Use slingshot for trajectory estimates / pseudotime
  - optimize cell type (Vinay!) estimation 
    - subsample "called" cell types, bootstrap and re-identify
    - ID "called" cell types with potential mis-calls
  - build prototype web app
  - 

![](https://github.com/davemcg/massive_integrated_eye_scRNA/blob/master/miescRNA.svg)
