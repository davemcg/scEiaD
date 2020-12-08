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
integrated_obj@meta.data$cluster <- meta[,2] %>% pull() %>% as.factor()
colnames(meta)[c(2,3)] <- c('cluster','subcluster')
umap <- umap %>% select(-cluster, -subcluster) %>% left_join(meta, by = 'Barcode')
seurat <- integrated_obj[, umap$Barcode]

if ('CellType_predict' %in% colnames(umap)){
	umap$CT <- umap$CellType_predict
} else {umap$CT <- umap$CellType}


quick_label_cluster <- function(umap){
	quick_label <-	umap %>%
	  #mutate(CT = gsub('Rod Bipolar Cells', 'Bipolar Cells', CT)) %>%
	  #mutate(CT = CT_predict) %>%
 	 group_by(cluster, CT) %>%
 	 filter(!is.na(CT), !is.na(cluster)) %>%
 	 summarise(Count = n(), x = mean(UMAP_1), y = mean(UMAP_2),
            Organism = list(unique(organism)),
            study_accession = list(unique(study_accession))) %>%
	  mutate(freq = Count / sum(Count)) %>%
 	 filter(freq > 0.25) %>%
  	ungroup() %>%
 	 group_by(cluster) %>%
  	top_n(3, -freq) %>%
 	 #filter(Count > 100) %>%
  	summarise(CT = paste0(CT, collapse = ' '),
            x = mean(x), y = mean(y),
            Count = sum(Count)) %>%
	filter(grepl('Precu|RPC|Gangli|Amacrine|Bipolar|Rods|Cones|Muller|Neurog|Horizon', CT)) %>%
	unique() %>% 
	ungroup() %>% 
	rowwise() %>% 
	mutate(seurat_cluster_CT = paste0(cluster, ': ', CT, collapse = ' ')) 	
	quick_label	
}

cut_down_objs <- function(org = 'all'){
	#umapRetina <- umap %>% filter(!is.na(CT)) %>% filter(grepl('Precu|RPC|Ganglia|Amacrine|Bipolar|Rods|Cones|Muller|Neurog|Horizon', CT))
	cl_quick <- quick_label_cluster(umap)
	umapRetinaCluster <- umap %>% filter(cluster %in% cl_quick$cluster) %>% left_join(cl_quick, by = 'cluster')
	if ('CellType_predict' %in% colnames(umap)){
		umapRetinaCluster <- umap %>% filter(CellType_predict %in% c('AC/HC_Precurs','Amacrine Cells','Astrocytes','Bipolar Cells','Cones','Early RPCs','Horizontal Cells','Late RPCs','Muller Glia','Neurogenic Cells','Pericytes','Photoreceptor Precursors','Retinal Ganglion Cells','Rod Bipolar Cells','Rods','RPCs') | is.na(CellType_predict)) %>% left_join(cl_quick, by = 'cluster')
	}
	if (org == 'all'){
		print('Using all cells')
		umapRetinaCluster = umapRetinaCluster
	} else {
		print(paste0('Using ', org, ' cells'))
		umapRetinaCluster = umapRetinaCluster %>% filter(organism == org)
	}
	sCT_CL = seurat[, umapRetinaCluster$Barcode]
	sCT_CL[["UMAP"]] <- CreateDimReducObject(embeddings =
                                umapRetinaCluster[,2:3] %>% as.matrix(),
                                key = "UMAP_",
                                assay = DefaultAssay(sCT_CL))
	list(umap = umapRetinaCluster, seurat = sCT_CL)
	
}	

run_sling <- function(seurat, group, reduction = 'scVI', ncell = 50000, start = NULL){
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
	} else {start <- grep(paste0('^', start), unique(colLabels(sce)), value = TRUE)
    }
	print(start)
	 ends <-  colLabels(sce) %>% 
					unique() %>% 
					enframe() %>% 
					filter(!grepl('RPC|Prec|Neuro', value) | grepl('Retinal Ganglio', value)) %>% 
					pull(value)
	print(ends)
	set.seed(90645)
	sceL <- sce[,sample(1:ncol(sce), ncell)]
	tic()
    sling  <- slingshot(sceL, 
						clusterLabels=colLabels(sceL), 
						reducedDim=toupper(reduction), 
						approx = 200, 
						start.clus = start,
						end.clus = ends) 
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

obj_cut <- cut_down_objs(org)
if (length(args) == 7) {
	start_clus <- args[7]
	sling <- run_sling(obj_cut$seurat, obj_cut$umap$seurat_cluster_CT, reduction, ncell = nrow(obj_cut$umap), start = start_clus)
} else {
	sling <- run_sling(obj_cut$seurat, obj_cut$umap$seurat_cluster_CT, reduction, ncell = nrow(obj_cut$umap))
}
umap_cut <- obj_cut$umap
save(umap, umap_cut, sling, file = out)
