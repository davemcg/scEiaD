# scEiaD

Goals: 
 - Harmonize all publicly **available** scRNA-seq datasets for eye (mostly back of the eye as of July 2020).
   - Also use Tabula Muris (pan mouse) data (July 2020)
 - Make {Shiny} app for data exploration
 - Stretch goals:
   - ID novel cell types
   - Species (in)consistent gene expression patterns across cell types
   - Splicing/transcript-level differences with well scRNA data
   - Pseudotime trajectories

   

Milestones:

 - Presented poster Genome Informatics 2019 (November). Managed to get a prototype integration with [fastMNN](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html).
   - Lior Pachter suggested using [scVI](https://scvi.readthedocs.io/en/stable/) to do the dimension reduction step
   - He also suggested using his/Sina's [kb python](https://github.com/pachterlab/kb_python) to run velocity estimates
 - Selected for talk at ARVO 2020 (May).
   - Talk cancelled due to coronavirus
 - "Completed" tasks:
   - Integration / Clustering:
     - Droplet only right now (see below for discussion)
     - Using CCA, combat, fastMNN, harmony, liger, magic, scanorama, scVI, insct , bbknn
       - CCA dropped (2020-07) for full integration as it runs too slowly / requires too much memory
       - still trying CCA with the well only data
     - Benchmarking with ARI, LISI, Silhouette
     - scVI better than the rest
     - Tested wide range of parameters (features, dimensions, knn, normalization) to pick best params for scVI
     - As of August 2020: 5000 features, 8 dims, knn 7
     - Hand compared walktrap with jaccard/louvain; the latter works better (visual inspection of cluster <-> umap)
     - Testing "subclustering" (higher resolution?) by running louvain again *within* a cluster
   - Diff testing:
     - Have dropped the Monocle based testing (2020-07), as it is too damn slow with ~1e6 cells and I can't get it to weigh the batches 
     - Have swapped over to pseudobulk based approach (https://osca.bioconductor.org/multi-sample-comparisons.html) with edgeR
    - Vinay/xgboost based CellType prediction (done 2020-07)
      - subsample "called" cell types, bootstrap and re-identify
      - ID "called" cell types with potential mis-calls
   - Droplet vs Well:
     - Originally I was trying to fully merge all data at once, but the well-based data was very problematic. I was going crazy. 
     - Then I read [Sina's](https://www.biorxiv.org/content/10.1101/2020.03.05.977991v3) mouse primary motor cortex paper. They have SMART-seq, 10X, (and MERFISH) data. Sina *independently* used each tech type and leveraged the largest advantage of SMART-seq: full(er) read length to get isoform specificity.
     - 2020-07 approach: three integrations:
       - onlyDroplet (major viz)
       - onlyWell (alt. viz)
       - universe (droplet + well) just to project celltypes onto well data
    - Internal (NEI) deploy of scAnthology Shiny app, April 2020
      - Thanks (mostly) to [scattermore](https://github.com/exaexa/scattermore)
      - [INSERT GIF?]
     - Cell Cycle ID (Done 2020-07. Stupid easy with Seurat.)
     - Doublet ID (done with scrublet and scran as of ~2020/04/25)
       - Have tried super hard to get solo to work, but it is too damn slow, which is surprising as it uses scVI
    - Tabula Muris
      - Integrated Tabula Muris Mouse Atlas data to provide non-retina cell type "comparison" (droplet only as of 2020-07)

 
