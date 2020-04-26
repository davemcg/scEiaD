# massive_integrated_eye_scRNA

Goal: Harmonize all publicly **available** scRNA-seq datasets for eye (mostly retina as of February 2020).

Milestones: 

<<<<<<< HEAD
 - Presented poster Genome Informatics 2019 (November). Managed to get a prototype integration with [fastMNN](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html).
   - Lior Pachter suggested using [scVI](https://scvi.readthedocs.io/en/stable/) to do the dimension reduction step
   - He also suggested using his/Sina's [kb python](https://github.com/pachterlab/kb_python) to run velocity estimates
 - Selected for talk at ARVO 2020 (May).
   - Talk cancelled due to coronavirus
 - "Completed" tasks:
   - Integration / Clustering:
     - Droplet only right now (see below for discussion)
     - Used CCA, combat, fastMNN, harmony, liger, magic, scanorama, scVI 
     - Benchmarking with ARI, LISI, Silhouette
     - scVI better than the rest
     - Tested wide range of parameters (features, dimensions, knn, normalization) to pick best params for scVI
     - As of April 2020: 2000 features, 200 dims, knn 7
     - Hand compared walktrap with jaccard/louvain; the latter works better (visual inspection of cluster <-> umap)
     - Testing "subclustering" (higher resolution?) by running louvain again *within* a cluster
   - Diff testing:
     - Monocle, two models
      - ~Cluster + percent.mt
      - ~Cluster + percent.mt + organism
   - Droplet vs Well:
     - Originally I was trying to fully merge all data at once, but the well-based data was very problematic. I was going crazy. 
     - Then I read [Sina's](https://www.biorxiv.org/content/10.1101/2020.03.05.977991v3) mouse primary motor cortex paper. They have SMART-seq, 10X, (and MERFISH) data. Sina *independently* used each tech type and leveraged the largest advantage of SMART-seq: full(er) read length to get isoform specificity.
     - So I'm going to crib off of this and run well/droplet separately (for integration) then post-hoc merge.
  - Internal (NEI) deploy of scAnthology Shiny app, April 2020
   - Thanks (mostly) to [scattermore](https://github.com/exaexa/scattermore)
   - [INSERT GIF?]
  - To do:
    - Velocity (see above)
    - optimize cell type (Vinay!) estimation
      - subsample "called" cell types, bootstrap and re-identify
      - ID "called" cell types with potential mis-calls
      - afterwords...can try to ID rarer subtypes???? 
      - check the Sanes more detailed RGC annotation to see if "subclustering" (above) resolves the **many** Sanes ID'ed RGC cell types
    - Finalize well-based processing (ARI/Silhouette/LISI calling is semi-borked on this smaller dataset)
    - "Integrate" well data into droplet ... somehow
    - Do isoform-level tx calling with well data
    - BIG ONE:
      - Process Mouse Atlas data to provide non-retina cell type "comparison" 
    - Continue to build out Shiny app
 
=======
 - Presented poster Genome Informatics 2019 (November). Managed to get a prototype integration with [fastMNN](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html). 
  - Lior Pachter suggested using [scVI](https://scvi.readthedocs.io/en/stable/) to do the dimension reduction step
  - He also suggested using his/Sina's [kb python](https://github.com/pachterlab/kb_python) to run velocity estimates
- As of Feb 2020 I've swapped over to scVI for the dimensionality reduction and tested a **wide** range of parameters and more rigorously assessed integration performance with ARI, LISI, and some home-spun methods as it was too onerous assessing the UMAP projections by eye. 
- I'm presenting this work at ARVO 2020 (May)! To do stuffs:
  - Use slingshot for trajectory estimates / pseudotime
  - optimize cell type (Vinay!) estimation 
    - subsample "called" cell types, bootstrap and re-identify
    - ID "called" cell types with potential mis-calls
    - afterwords...can try to ID rarer subtypes???? perhaps ID groups/clusters/sections with "fuzzy" confidence in annotation?
    - check the sanes more detailed RGC annotation
  - build prototype Shiny web app
     - quick viz (scattermore)
     - summary info (heatmap?) by cell type / species / time
     - give info on # of indepedent studies/batches supporting cluster / group
    

![](https://github.com/davemcg/massive_integrated_eye_scRNA/blob/master/miescRNA.svg)
>>>>>>> parent of 2509e55... Update README.md
