#Use >= R/4.0
library(tidyverse)
library(scran)
library(slingshot)
library(Seurat)
library(tictoc)

args = commandArgs(trailingOnly=TRUE)
load(args[1]) # load('cluster/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__knn0.6.cluster.Rdata') 
load(args[2]) # load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.001__nneighbors500.umapFilter.predictions.Rdata')
load(args[3]) # load('seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__counts__TabulaDroplet__batch__scVI__dims8__preFilter.seuratV3.Rdata')
method = args[4]
org = gsub('_', ' ', args[5])
out <- args[6]
subdivide_CT_by_clusters <- args[7]
integrated_obj@meta.data$cluster <- meta[,2] %>% pull() %>% as.factor()
colnames(meta)[c(2,3)] <- c('cluster','subcluster')
umap <- umap %>% select(-cluster, -subcluster) %>% left_join(meta, by = 'Barcode')
seurat <- integrated_obj[, umap$Barcode]

print("ARGS")
print(args)
print("")

if ('CellType_predict' %in% colnames(umap)){
	umap$CT <- umap$CellType_predict
} else {umap$CT <- umap$CellType}


cut_down_objs <- function(umap, org = 'Homo sapiens'){
  umapRetina = umap %>% filter(organism == org,
                  is.na(TabulaMurisCellType_predict),
                  CellType_predict %in% c('AC/HC_Precurs','Amacrine Cells','Astrocytes',
                                           'Bipolar Cells','Cones','Early RPCs','Horizontal Cells', 'RPE',
                                           'Late RPCs','Muller Glia','Neurogenic Cells','Photoreceptor Precursors',
                                           'Retinal Ganglion Cells','Rod Bipolar Cells','Rods','RPCs') | is.na(CellType_predict))
  
  sCT_CL = seurat[, umapRetina$Barcode]
  sCT_CL[["UMAP"]] <- CreateDimReducObject(embeddings =
                                             umapRetina[,2:3] %>% as.matrix(),
                                           key = "UMAP_",
                                           assay = DefaultAssay(sCT_CL))
  
  list(umap = umapRetina, seurat = sCT_CL)
  
}

run_sling <- function(seurat, group, reduction = 'scVI', ncell = 50000, start = NULL, omega = Inf, sling_guess_ends = FALSE){
  sce <- as.SingleCellExperiment(seurat)
  sce$group <- group
  colLabels(sce) <- sce$group
  if (length(umap$organism %>% unique()) == 1 & umap$organism %>% unique()  == 'Mus musculus'){
    grep_against <- 'bloop'
  } else {grep_against <- '^RPC'
  }
  if (is.null(start)){
    start <- colLabels(sce) %>%
      table() %>%
      enframe() %>%
      arrange(-value) %>%
      filter(!grepl(grep_against, name)) %>%
      filter(grepl('RPC', name)) %>%
      head(1) %>%
      pull(name)
  } 
  #else {start <- grep(paste0('^', start, ':'), unique(colLabels(sce)), value = TRUE)
  #}
  print("Start Cluster")
  print(start)
  print('')

  ends <-  colLabels(sce) %>%
    unique() %>%
    enframe() %>%
    filter(!grepl('RPC|Prec|Neuro', value) ) %>%
    pull(value)
  
  intermediate <- colLabels(sce) %>%
    unique() %>%
    enframe() %>%
    filter(!value %in% ends) %>% 
    pull(value)
  
  set.seed(90645)
  sceL <- sce[,sample(1:ncol(sce), ncell)]
  tic()
  if (sling_guess_ends == FALSE){
    print("intermediate groups")
    print(intermediate %>% sort())
    print("End clusters")
    print(ends %>% sort())
    print('')
    sling  <- slingshot(sceL,
                        clusterLabels=colLabels(sceL),
                        reducedDim=toupper(reduction),
                        approx = 200,
                        omega = omega,
                        start.clus = start,
                        end.clus = ends)
  } else {
    print('Sling will guess ends!!')
    sling  <- slingshot(sceL,
                        clusterLabels=colLabels(sceL),
                        reducedDim=toupper(reduction),
                        approx = 200,
                        omega = omega,
                        start.clus = start)
  }
  toc()
  lineage <- getLineages(reducedDim(sceL, type = 'UMAP'),  colLabels(sceL),  start.clus = start, end.clus = ends)
  out <- list()
  out$sling <- sling
  out$lineage <- lineage
  out$start <- start
  out$ends <- ends
  out$embedded <- embedCurves(sling, "UMAP")
  out
}

if (method == 'CCA'){
  reduction <- 'pca'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
} else if (method == 'none'){
  reduction <- 'pca'
} else if (method == 'combat'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
} else if (method == 'liger'){
  reduction <- 'iNMF'
} else if (grepl('scVI', method)){
   reduction <- 'scVI'
} else {
  print("GUESSING!")
  reduction <- method
}

rm(integrated_obj)

obj_cut <- cut_down_objs(umap, org = 'Homo sapiens')
# now's the time to pick a start cluster
#quick_label_cluster(obj_cut$umap) %>% data.frame()

if (length(args) == 8) {
	print('Start  cluster given as user input')
	start_clus <- args[8]
	sling <- run_sling(obj_cut$seurat, obj_cut$umap$seurat_cluster_CT, reduction, ncell = nrow(obj_cut$umap), start = start_clus)
} else {
	sling <- run_sling(obj_cut$seurat, obj_cut$umap$seurat_cluster_CT, reduction, ncell = nrow(obj_cut$umap))
}

umap_cut <- obj_cut$umap



#### diff testing time
library(tidyverse)
library(scran)
library(TSCAN)
library(slingshot)
library(tictoc)
args = commandArgs(trailingOnly=TRUE)
#load(args[1]) # slingshot obj

pt_tester <- function(sce, pseudoTime, block){
  pT <- colData(sce)[,pseudoTime]
  block <- colData(sce)[,block]
  # remove NA pseudotimes
  pT_clean <- pT[!is.na(pT)]
  block_clean <- block[!is.na(pT)]
  # aggregate low n (<100) blocks into 1
  low_n_blocks <- block_clean %>% table() %>% enframe() %>% filter(value < 200) %>% pull(name)
  new_name <- paste(low_n_blocks, collapse = '__')
  block_clean[block_clean %in% low_n_blocks] <- new_name
  # diff test
  diff_PT <- TSCAN::testPseudotime(sce[,!is.na(pT)],
                                   pseudotime = pT_clean,
                                   block = block_clean,
                                   BPPARAM = BiocParallel::MulticoreParam(12))
  diff_PT
}

sce <- sling$sling
diffPT <- list()
for (i in colnames(colData(sce)) %>% grep('slingP', .,  value = TRUE)){
  print(i)
  tic()
  diffPT[[i]] <- try({pt_tester(sce, i, 'study_accession')})
  toc()
}

save(umap, umap_cut, sling, diffPT, file = out)
