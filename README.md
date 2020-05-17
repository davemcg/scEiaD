# scAnthology

Goals: 
 - Harmonize all publicly **available** scRNA-seq datasets for eye (mostly retina as of February 2020).
 - Make {Shiny} app for data exploration
 - Stretch goals:
   - ID novel cell types
   - Species (in)consistent gene expression patterns across cell types
   - Splicing/transcript-level differences with well scRNA data

   

Milestones:

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
     - Monocle, four models
      - ~Cluster + percent.mt + batch
      - ~Cluster + percent.mt + batch + organism
      - ~CellType_predict + percent.mt + batch
      - ~CellType_predict + percent.mt + batch + organism    
   - Droplet vs Well:
     - Originally I was trying to fully merge all data at once, but the well-based data was very problematic. I was going crazy. 
     - Then I read [Sina's](https://www.biorxiv.org/content/10.1101/2020.03.05.977991v3) mouse primary motor cortex paper. They have SMART-seq, 10X, (and MERFISH) data. Sina *independently* used each tech type and leveraged the largest advantage of SMART-seq: full(er) read length to get isoform specificity.
     - So I'm going to crib off of this and run well/droplet separately (for integration) then post-hoc merge.
  - Internal (NEI) deploy of scAnthology Shiny app, April 2020
    - Thanks (mostly) to [scattermore](https://github.com/exaexa/scattermore)
    - [INSERT GIF?]
  - To do:
    - [Velocity](https://bustools.github.io/BUS_notebooks_R/velocity.html) 
    - Cell Cycle ID
    - Doublet ID (done with scrublet as of ~2020/04/25)
    - Diff testing
      - Diff testing on subclustering
      - Implement [pseudobulk](https://osca.bioconductor.org/multi-sample-comparisons.html#differential-expression-between-conditions) by cell type / cluster?
      - Run [Wilcox test](https://osca.bioconductor.org/marker-detection.html#using-the-wilcoxon-rank-sum-test) to get cluster / etc AUC and combine with Monocle testing in unified DF (done with scran as of ~2020/04/25)
    - optimize cell type (Vinay!) estimation
      - subsample "called" cell types, bootstrap and re-identify
      - ID "called" cell types with potential mis-calls
      - afterwords...can try to ID rarer subtypes???? 
      - check the Sanes more detailed RGC annotation to see if "subclustering" (above) resolves the **many** Sanes ID'ed RGC cell types
    - Finalize well-based processing (ARI/Silhouette/LISI calling is semi-borked on this smaller dataset)
    - "Integrate" well data into droplet ... somehow
      - easiest is to not use for viz, but just add to celltype-based tables
      - cluster2pattern from [projectR](https://www.bioconductor.org/packages/release/bioc/vignettes/projectR/inst/doc/projectR.pdf)?
    - Do isoform-level tx calling with well data
      - Crib from [Sina's](https://www.biorxiv.org/content/10.1101/2020.03.05.977991v3) paper
    - BIG ONE:
      - Process Mouse Atlas data to provide non-retina cell type "comparison" 
    - Continue to build out Shiny app
 
