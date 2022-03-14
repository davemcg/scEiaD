library(glue)
library(yaml)

make_seurat_obj <- function(m,
                            cell_info,
                            split.by = 'study_accession',
                            nfeatures = n_features,
                            keep_well = TRUE,
                            keep_droplet = TRUE,
                            mito_geneids,
							lengthCor = FALSE, 
							dont_use_well_for_FVF = FALSE,
							only_use_human_for_FVF = FALSE,
							droplet_platform = c('DropSeq', '10xv2', '10xv3'),
							well_platform = c('SMARTSeq_v2', 'SCRBSeq', 'SMARTerSeq_v3', 'SMARTSeq_v4'),
							HVG = NA,
							nFeature_RNA_cutoff = 400,
							custom_nFeature = TRUE
                            ){
  nFeature_RNA_cutoff <- 450
  print(paste('nFeature_RNA_cutoff is', nFeature_RNA_cutoff))
  if (keep_well){
      well_m <- m[,cell_info %>% filter(value %in% colnames(m), Platform %in% well_platform) %>% pull(value)]
     seurat_well <- CreateSeuratObject(well_m)
  }
  if (keep_droplet & custom_nFeature){
     # custom handling for wacky E12 chick sanes data
	 # trying a stricter cutoff as way too many cells (getting ~200k and their paper used ~50 ) are being kept
     # SRP362101 getting nearly 200k and they only should have around 15

	droplet_subset_maker <- function(m, samples, cutoff){
	 	 droplet_sub <- m[,cell_info %>% filter(sample_accession %in% samples, value %in% colnames(m)) %>% pull(value)]
		 seurat_sub <- CreateSeuratObject(droplet_sub)
		 seurat_sub@meta.data %>% as_tibble(rownames = 'Barcode') %>% filter(nFeature_RNA < cutoff) %>% pull(Barcode)
	}	 
	 
	 e12_chick_bc_rm <- droplet_subset_maker(m,  c('SRS7483328','SRS7483329','SRS7483330','SRS7483331'), 1000)
     SRP362101_bc_rm <- droplet_subset_maker(m, samples = c("SRX14336032","SRX14336031","SRX14336030","SRX14336029","SRX14336028","SRX14336027","SRX14336026","SRX14336025","SRX14336024","SRX14336023","SRX14336022","SRX14336021","SRX14336020","SRX14336019","SRX14336018","SRX14336017","SRX14336016","SRX14336015","SRX14336014","SRX14336013","SRX14336012","SRX14336011","SRX14336010","SRX14336009","SRX14336008","SRX14336007","SRX14336006","SRX14336005","SRX14336004","SRX14336003","SRX14336002","SRX14336001","SRX14336000","SRX14335999","SRX14335998","SRX14335997","SRX14335996","SRX14335995","SRX14335994","SRX14335993","SRX14335992","SRX14335991","SRX14335990","SRX14335989","SRX14335988","SRX14335987","SRX14335986","SRX14335985","SRX14335984","SRX14335983","SRX14335982","SRX14335981","SRX14335980","SRX14335979","SRX14335978","SRX14335977","SRX14335976","SRX14335975","SRX14335974","SRX14335973","SRX14335972","SRX14335971","SRX14335970","SRX14335969","SRX14335968","SRX14335967","SRX14335966","SRX14335965","SRX14335964","SRX14335963","SRX14335962","SRX14335961","SRX14335960","SRX14335959","SRX14335958","SRX14335957","SRX14335956","SRX14335955","SRX14335954","SRX14335953","SRX14335952","SRX14335951","SRX14335950","SRX14335949","SRX14335948","SRX14335947","SRX14335946","SRX14335945","SRX14335944","SRX14335943","SRX14335942","SRX14335941","SRX14335940","SRX14335939","SRX14335938","SRX14335937","SRX14335936","SRX14335935","SRX14335934","SRX14335933","SRX14335932","SRX14335931","SRX14335930","SRX14335929","SRX14335928","SRX14335927","SRX14335926","SRX14335925","SRX14335924","SRX14335923","SRX14335922","SRX14335921","SRX14335920","SRX14335919","SRX14335918","SRX14335917","SRX14335916","SRX14335915","SRX14335914","SRX14335913","SRX14335912","SRX14335911","SRX14335910","SRX14335909","SRX14335908","SRX14335907","SRX14335906","SRX14335905","SRX14335904","SRX14335903","SRX14335902","SRX14335901","SRX14335900","SRX14335899","SRX14335898","SRX14335897","SRX14335896","SRX14335895","SRX14335894","SRX14335893","SRX14335892","SRX14335891","SRX14335890","SRX14335889"), cutoff = 800)

	 #  use seurat_e12_chick and seurat_SRP362101 as cells to skip in this line:
     droplet_m <- m[,cell_info %>% filter(value %in% colnames(m), !value %in% c(e12_chick_bc_rm, SRP362101_bc_rm), Platform %in% droplet_platform) %>% pull(value)]
     
	 seurat_droplet <- CreateSeuratObject(droplet_m)
	 seurat_droplet <-  subset(seurat_droplet, subset = nFeature_RNA > 450)
     
	 print('seurat droplet obj made')
	 #droplet_hs <- m[,cell_info %>% filter(value %in% colnames(m), organism %in% c('Homo sapiens')) %>% pull(value)]
	 #hs_seurat_droplet <- CreateSeuratObject(droplet_hs)
  }

  if (keep_droplet & !custom_nFeature){
    droplet_m <- m[,cell_info %>% filter(value %in% colnames(m), Platform %in% droplet_platform) %>% pull(value)] 
	seurat_droplet <- CreateSeuratObject(droplet_m)
	 seurat_droplet <-  subset(seurat_droplet, subset = nFeature_RNA > 450)
     
	 print('seurat droplet obj made')
  }
  # FILTER STEP!!!!
  # keep cells with < 10% mito genes, and more than 200 and less than 3000 detected genes for UMI
  # for well, drop the 3000 gene top end filter as there shouldn't be any droplets
  if (keep_well & !lengthCor){
    print('No Length Correction')
    seurat_well <- subset(seurat_well, subset = nFeature_RNA > 500)
  } else if (keep_well && lengthCor) {
	  mm <- FALSE
	  if ('mouse' %in% set){
		mm <- TRUE
	  }
      source(glue('{git_dir}/src/extract_gene_length.R'))
    	geneL_mm <- gene_length( 'references/gtf/mm-mus_musculus_anno.gtf.gz', mm)
    	geneL_hs <-  gene_length('references/gtf/hs-homo_sapiens_anno.gtf.gz')
    	well_hs <- cell_info %>% filter(organism == 'Homo sapiens') %>% pull(value)
        well_mm <- cell_info %>% filter(organism == 'Mus musculus') %>% pull(value)
    	mat <- seurat_well@assays$RNA@counts
    	hs_mat <- mat[, well_hs[well_hs %in% colnames(mat)]]
    	mm_mat <-  mat[, well_mm[well_mm %in% colnames(mat)]]
    	hs_mat_cor <- hs_mat / geneL_hs[row.names(hs_mat)] * median(geneL_hs)
    	mm_mat_cor <- mm_mat / geneL_mm[row.names(mm_mat)] * median(geneL_mm)
		hs_mm <- cbind(hs_mat_cor, mm_mat_cor)
    	mat_cor <- hs_mm[, colnames(hs_mm) %in% colnames(mat)] %>% round()
		mat_cor[is.na(mat_cor)] <- 0
  	  seurat_well <- CreateSeuratObject(mat_cor)
      seurat_well <- subset(seurat_well, subset = nFeature_RNA > 500)
      print('Length Correction for Well DONE!')
  }
  # cells to keep
  if (!keep_well & keep_droplet){
    print('removing well')
    cells_to_keep <- row.names(seurat_droplet@meta.data)
  } else if (keep_well & !keep_droplet) {
    print('removing droplet')
    cells_to_keep <- row.names(seurat_well@meta.data)
  } else {
    print('keeping well and droplet')
    cells_to_keep <- c(row.names(seurat_droplet@meta.data), row.names(seurat_well@meta.data))
  }


  m_filter <- m[,cells_to_keep]


  seurat_m <- CreateSeuratObject(m_filter)
  mito_geneids_present <- mito_geneids[mito_geneids%in% rownames(seurat_m)]
  seurat_m[["percent.mt"]] <- PercentageFeatureSet(seurat_m, features = mito_geneids_present)
  seurat_m@meta.data$batch <- left_join(seurat_m@meta.data %>%
                                          row.names() %>% enframe(),
                                        cell_info, by = 'value') %>%
    pull(batch)
  seurat_m@meta.data$study_accession <- left_join(seurat_m@meta.data %>%
                                                    row.names() %>% enframe(),
                                                  cell_info, by = 'value') %>%
    pull(study_accession)
  seurat_m@meta.data$Age <- left_join(seurat_m@meta.data %>%
                                        row.names() %>% enframe(),
                                      cell_info, by = 'value') %>%
    pull(Age)
  seurat_m@meta.data$TechType <- left_join(seurat_m@meta.data %>%
                                        row.names() %>% enframe(),
                                      cell_info, by = 'value') %>%
    pull(Platform)

  ## account for any batches with less than 3 samples 
  batch_count <- table(seurat_m$batch) 
  bad_batches <- names(batch_count[batch_count<3])
  if(length(bad_batches) > 0){
    print('Removing the following batches which had ncells <3')
    print(bad_batches)
    good_cells <- filter(seurat_m@meta.data, !batch %in% bad_batches ) %>% rownames
    seurat_m <- seurat_m[,good_cells]
  }
 	
  # MT filtering
  # 10xv3 has approx 2x the mito counts of 10xv2, so giving that platform a diff cutoff
  # https://kb.10xgenomics.com/hc/en-us/articles/360026501692-Do-we-see-a-difference-in-expression-profile-of-3-Single-Cell-v3-chemistry-as-compared-to-v2-chemistry
  mt_cells_keep <- (seurat_m@meta.data$percent.mt < 10 & seurat_m@meta.data$TechType != '10xv3') | (seurat_m@meta.data$percent.mt < 20 & seurat_m@meta.data$TechType == '10xv3')
  seurat_m <- seurat_m[, mt_cells_keep]
 
  # scale data and regress
  seurat_m <- NormalizeData(seurat_m)
  # find var features
  seurat_m <- FindVariableFeatures(seurat_m, nfeatures = nfeatures, selection.method = 'vst')
  if (dont_use_well_for_FVF == TRUE) {
	seurat_mDrop <- FindVariableFeatures(seurat_droplet, nfeatures = nfeatures, selection.method = 'vst')
    VariableFeatures(seurat_m) <- VariableFeatures(seurat_mDrop)
  }

  if (only_use_human_for_FVF == TRUE) {
	print('only use human for FVF')
    droplet_hs <- m[,cell_info %>% filter(value %in% colnames(seurat_m), Platform %in% droplet_platform, organism %in% c('Homo sapiens')) %>% pull(value)]
    hs_seurat_droplet <- CreateSeuratObject(droplet_hs)
	hs_seurat_droplet  <- FindVariableFeatures(hs_seurat_droplet, nfeatures = nfeatures, selection.method = 'vst')
    VariableFeatures(seurat_m) <- VariableFeatures(hs_seurat_droplet)
  }  
  # don't use mito genes
  var_genes <- grep('^MT-', seurat_m@assays$RNA@var.features, value = TRUE, invert = TRUE)
 
  if (!is.na(HVG)) {
		# ensure the HVG are actually present
		HVGthere = HVG[HVG %in% row.names(seurat_m)]
		VariableFeatures(seurat_m) <- HVGthere[!is.na(HVGthere[1:nfeatures])]
  } 

  if (transform == 'standard'){
    print(paste0('Running lib.size and log correction, splitting by ', split.by))

    data <- seurat_m@assays$RNA@counts %>% t()
    data <- data[,var_genes]
    library_size <- Matrix::rowSums(data)
    median_transcript_count <- stats::median(library_size)
    data_norm <- median_transcript_count * data / library_size
    data_norm <- t(data_norm)
    data_norm <- log(data_norm + 1)
    seurat_m@assays$RNA@scale.data <- data_norm %>% as.matrix()
     #seurat_m <- ScaleData(seurat_m,
      #                    features = var_genes,
       #                   split.by = split.by,
        #                  do.center = TRUE,
         #                 do.scale = TRUE,
          #                verbose = TRUE,
           #               vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

    seurat_m <- RunPCA(seurat_m, npcs = 100)
  }
  seurat_m
}


# build SCT based seurat obj
seurat_sct <- function(seurat_list){
  ## tryCatch for SCT
  trySCTransform <- function(x){
    tryCatch(
      expr = {
        message("Successful SCTransform")
        return(SCTransform(x, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")))
      },
      error = function(e){
        message("Failed SCTransform")
        print(e)
        return(NULL)
      }
    )
  }
  for (i in names(seurat_list)){
    DefaultAssay(seurat_list[[i]]) <- 'RNA'
    seurat_list[[i]] <- trySCTransform(seurat_list[[i]])
  }

  # remove sets with less than 500 cells, which will somehow(?) destory SCT - based integration performance
  low_n <- c()
  for (i in names(seurat_list)){
    if (ncol(seurat_list[[i]]) < 500){
      low_n <- c(low_n, i)
    }
  }
  if (length(low_n) > 0){
    seurat_list[low_n] <- NULL
  }

  study_data_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = n_features, verbose = FALSE)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = study_data_features, verbose = FALSE)

  # have to apply PCA after SelectIntegrationFeatures and PrepSCTIntegration
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- RunPCA(x, features = study_data_features, verbose = FALSE)
  })
  list(seurat_list = seurat_list, study_data_features = study_data_features)
}

#cran normalization
scran_norm <- function(seurat_obj = seurat__standard, split.by = 'batch'){
  var_features <- grep('^MT-', seurat_obj@assays$RNA@var.features, invert =TRUE, value = TRUE)
  seurat_list <- SplitObject(seurat_obj, split.by = covariate)
  print('Beginning scran norm')
  # list of seurat objects
  sce_list <- list()
  for (obj in names(seurat_list)){
    print(obj)
    print(seurat_list[[obj]] %>% dim())
    if (obj == 'SRP131661_10xv2_3-M-8'){
	  seurat_list[[obj]] <- NULL
	} else if (ncol(seurat_list[[obj]]) > 100){
      sce_list[[obj]] <- SingleCellExperiment(assays = list(counts = as.matrix(x = seurat_list[[obj]]$RNA@data)))
      clusters = quickCluster(sce_list[[obj]], min.size=100)
      sce_list[[obj]] = computeSumFactors(sce_list[[obj]], cluster=clusters)
      sce_list[[obj]] = scater::logNormCounts(sce_list[[obj]], log=TRUE)
      #sce_list[[obj]] = normalize(sce_list[[obj]], return_log = FALSE)
      print(summary(sizeFactors(sce_list[[obj]])))
      # seurat_list[[obj]]$RNA@data = as.sparse(log(x = assay(sce_list[[obj]], "normcounts") + 1))
      seurat_list[[obj]]$RNA@data = assay(sce_list[[obj]], "logcounts") %>% as.sparse()
    } else {seurat_list[[obj]] <- NULL} # remove obj with less than 50 cells
  }
  # merge back into one seurat obj
  merged <- merge(x = seurat_list[[1]], y = seurat_list[2:length(x = seurat_list)])
  merged$RNA@scale.data = merged$RNA@data[var_features,] %>% as.matrix()
  merged@assays$RNA@var.features <- var_features
  # re do PCA
  merged <- RunPCA(merged, npcs = 100, features = var_features)
  merged
}



#' Performs L1 normalization on input data such that the sum of expression
#' values for each cell sums to 1, then returns normalized matrix to the metric
#' space using median UMI count per cell effectively scaling all cells as if
#' they were sampled evenly.
#' @param seurat_object
#' @return seurat_obj
#' 2 dimensional array with normalized gene expression values
#' @import Matrix
#' @import dplyr
#'
#' @export
library.size.normalize <- function(seurat_obj, sqrt = FALSE, verbose=FALSE) {
  if (verbose) {
    message(paste0(
      "Normalizing library sizes for ",
      nrow(data), " cells"
    ))
  }
  vfeatures <- grep('^MT-', seurat_obj@assays$RNA@var.features, invert =TRUE, value = TRUE)
  data <- seurat_obj@assays$RNA@counts %>% t()
  data <- data[,vfeatures]
  library_size <- Matrix::rowSums(data)
  median_transcript_count <- stats::median(library_size)
  data_norm <- median_transcript_count * data / library_size
  data_norm <- t(data_norm)
  if (sqrt) {data_norm <- sqrt(data_norm)}
  seurat_obj@assays$RNA@scale.data <- data_norm %>% as.matrix()
  seurat_obj <- RunPCA(seurat_obj, features = vfeatures)
  seurat_obj
}


