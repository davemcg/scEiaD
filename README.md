# massive_integrated_eye_scRNA

Goal: identify all publicly available scRNA-seq datasets for eye (mostly retina right as of June 2019).

Approach:
  - search GEO for "retina single RNA" and "retina scRNA"
  - collect GEO IDs from NCBI search results (e.g. 200130636, 200119343, 200119274, 200112507)
  - grab some metadata by by creating NCBI search term
    - e.g. https://www.ncbi.nlm.nih.gov/gds/?term=200130636+200119343+200119274+200112507
  - get SRA project ID
  - use omicidx to get SRA accession nums
  - clean metadata
  - download with fastq-dump 
  - counts with bustools - kallisto
  - merge all with liger
    - https://www.biorxiv.org/content/10.1101/459891v1
  - analyze:
    - species
    - cell types
    - trajectory
