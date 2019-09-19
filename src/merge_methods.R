# run integration methods that support Seurat objects directly
args <- commandArgs(trailingOnly = TRUE)

method = args[1]
covariate = args[2]
load(args[3])

# crazy section to deal with that fact I have scanorama in a conda environment,
# but many of the seurat wrapped integration tools can't be installed in conda
# without crazy effort (e.g liger...as it needs to be compiled in C)
library(tidyverse)
library(Seurat)
if (method != 'scanorama'){
  library(SeuratWrappers)
  library(harmony)
  library(batchelor)
} else {
  library(reticulate)
  scanorama <- import('scanorama')
}

run_integration <- function(seurat_obj, method, covariate = 'study_accession'){
  # covariate MUST MATCH what was used in build_seurat_obj.R
  # otherwise weird-ness may happen
  # the scaling happens at this level
  # e.g. DO NOT use 'batch' in build_seurat_obj.R then 'study_accession' here
  if (method == 'CCA'){
    # identify batches wiht really low cell counts (<100) to exclude
    obj <- seurat_obj
	obj@meta.data$split_by <- obj@meta.data[,covariate]
    splits_to_remove <- obj@meta.data %>% 
      dplyr::group_by(split_by)%>% 
      summarise(Count = n())  %>% 
      filter(Count < 100) %>% 
      pull(split_by)
	if (length(splits_to_remove >= 1)) { 
    	obj <- subset(obj, subset = split_by %in% splits_to_remove, invert = TRUE)
    }
    seurat_list <- SplitObject(obj, split.by = covariate)
    anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20)
    obj <- IntegrateData(anchorset = anchors, verbose = TRUE)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, npcs = 100)
  } else if (method == 'fastMNN'){
    ## uses list of seurat objects (each obj your "covariate")
    seurat_list <- SplitObject(seurat_obj, split.by = covariate)
    obj <- RunFastMNN(object.list = seurat_list)
  } else if (method == 'harmony'){
    ## uses one seurat obj (give covariate in meta.data to group.by.vars)
    obj <- RunHarmony(seurat_obj, group.by.vars = covariate,
                      max.iter.harmony = 10, 
                      epsilon.harmony = -Inf)
  } else if (method == 'liger'){
    ## like harmony above, give one seurat obj
    ## NMF requires POSITIVE values
    ## so, UNLIKE, the tutorial set the lowest value to 0.1 in the matrix
    var_genes <- grep('^MT-', seurat_obj@assays$RNA@var.features, value = TRUE, invert = TRUE)
    obj <- seurat_obj
    # identify batches wiht really low cell coutns (<100) to exclude
    obj@meta.data$split_by <- obj@meta.data[,covariate]
    splits_to_remove <- obj@meta.data %>% 
      dplyr::group_by(split_by)%>% 
      summarise(Count = n())  %>% 
      filter(Count < 100) %>% 
      pull(split_by) 
	if (length(splits_to_remove >= 1)) { 
    	obj <- subset(obj, subset = split_by %in% splits_to_remove, invert = TRUE)
    }
    #obj <-  ScaleData(seurat_obj, split.by = covariate, vars.to.regress = var_genes, do.center = FALSE)
    obj@assays$RNA@scale.data <- obj@assays$RNA@scale.data - min(obj@assays$RNA@scale.data) + 0.1 # <- yeah this is hacky but I think OK...
    # ...the alternative would be to re-run from scratch with raw counts, which would mean 
    # liger would be getting differently scaled values
    # than the other methods...
    obj <- RunOptimizeALS(obj, 
                          k = 20, 
                          lambda = 5, 
                          split.by = covariate)
    obj <- RunQuantileAlignSNF(obj, split.by = "Method")
    # finally re-do scaledata in case I use it in the future...I prob won't be expecting
    # it to not be centered and with no neg values....
    obj <- ScaleData(obj,  
                     features = var_genes,
                     do.center = TRUE,
                     do.scale = TRUE,
                     vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
  } else if (method == 'scanorama'){
    # scanorama can return "integrated" and/or "corrected" data
    # authors say that the "integrated" data is a low-dimension (100) representation
    # of the integration, which is INTENDED FOR PCA/tSNE/UMAP!!!
    # the corrected data returns all of the sample x gene matrix with batch
    # corrected values
    assay <- 'RNA'
    d <- list() # d is lists of scaled Data
    g <- list() # g is lists of gene names for each matrix in d
    seurat_list <- SplitObject(seurat_obj, split.by = covariate)
    # build name-less list of expression values
    # and name-less list of gene names
    # if you give a list with names then reticulate / python integration explodes
    # well actully just the name are given to python...which took me
    # forever to figure out. stupid annoying. 
    for (i in seq(1, length(seurat_list))){
      # print(i);
      var_genes <- grep('^MT-', seurat_list[[i]]@assays$RNA@var.features, value = TRUE, invert = TRUE)
      d[[i]] <- t((seurat_list[[i]]@assays[[assay]]@scale.data[var_genes,])) %>% as.matrix(); 
      d[[i]][is.na(d[[i]])] <- 0; 
      g[[i]] <- colnames(d[[i]]) 
    }
    print('Run scanorama')
    integrated.corrected.data <- 
      scanorama$correct(d, g , return_dimred=TRUE, return_dense=TRUE)
    print('Scanorama done')
    # get the cell names back in as they aren't returned for some reason (reticulate?)
    for (i in seq(1, length(d))){
      # first is in the integrated data (dim reduced for UMAP, etc)
      row.names(integrated.corrected.data[[1]][[i]]) <- row.names(d[[i]])
    }
    
    # glue reduced matrix values into seurat for later UMAP, etc
    scanorama_mnn <- Reduce(rbind, integrated.corrected.data[[1]])
    colnames(scanorama_mnn) <- paste0("scanorama_", 1:ncol(scanorama_mnn))
    obj <- seurat_obj
    obj[["scanorama"]] <- CreateDimReducObject(embeddings = scanorama_mnn, key = "scanorama_", assay = DefaultAssay(obj))
    
  } else {
    print('Supply either CCA, fastMNN, harmony, liger, or scanorama as a method')
    NULL
  }
  obj
}

extract_umap_for_plotting <- function(integrated_obj){
  load('Mus_musculus_cell_info_labelled.Rdata')
  orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
  umap <- Embeddings(integrated_obj[['umap']]) %>% as_tibble(rownames = 'Barcode') %>% 
    left_join(., orig_meta) %>% 
    left_join(., cell_info_labels %>% select(value:Paper) %>% rename(Barcode = value))
}

integrated_obj <- run_integration(seurat__standard, method, covariate)
save(integrated_obj, file = args[4], compress = FALSE)
