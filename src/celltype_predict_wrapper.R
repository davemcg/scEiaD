library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

load(args[1]) # seurat obj post integration
load(args[2]) # umap metadata

conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
working_dir = Sys.getenv('SCIAD_WORKING_DIR')

out <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode') %>% left_join(umap, by = 'Barcode')
out_tm <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode') %>% left_join(umap, by = 'Barcode') %>% filter(study_accession == 'SRP131661')
out_tm$CellType <- out_tm$TabulaMurisCellType
remove_low <- out_tm %>% group_by(CellType) %>% count() %>% filter(n < 40) %>% pull(CellType)
out_tm <- out_tm %>% filter(!CellType %in% remove_low)
# fantom5 tf
#f5 <- read_tsv('~/git/massive_integrated_eye_scRNA/data/fantom5_tf.tsv')
#HVTF <- f5$Symbol[(f5$Symbol %in% integrated_obj@assays$RNA@var.features)]
#counts <- integrated_obj@assays$RNA@counts[HVTF,]
#cpm <- RelativeCounts(counts, scale.factor= 1e6) %>% as.matrix() %>% t()

#out <- cbind(out, cpm)

# remove celltype outliers
rm_outlier <- function(out, TM = FALSE){
	bc_retain <- list()
	for (i in out$CellType %>% unique){
		print(i)
		temp <- out %>% filter(CellType == i) %>% select(contains('scVI_')) %>% as.matrix()
		row.names(temp) <- out %>% filter(CellType == i) %>% pull(Barcode)
		scVI_mean <- colMeans(temp)
		euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
		D <- apply(as.matrix(temp), 1, function(x) euc_dist(scVI_mean, x))
		cutoff = sd(D) * 4
		bc_retain[[i]] <- D[D < cutoff] %>% names()
	}
	out %>% filter(Barcode %in% (bc_retain %>% unlist()))
}

full_out <- out
full_tm <- out_tm
out <- rm_outlier(out)
out_tm <- rm_outlier(out_tm)

run_xgboost_py <- function(embeddings, full_embeddings, temp_name, tm = FALSE, probThresh = 0.8){
	rand_num <- sample(1e5:1e6, 1)
	embeddings_file <- paste0(rand_num, '_', temp_name)
	full_embeddings_file <- paste0(rand_num, '_FULL', temp_name)
	write_tsv(embeddings, path = embeddings_file)
	write_tsv(full_embeddings, path = full_embeddings_file)
	write_features_file <- paste0(rand_num, '_features.txt')
	if (!tm) {
		features <- c(grep('scVI', colnames(embeddings), value = TRUE),
						'UMAP_1', 'UMAP_2', 'nCount_RNA', 'nFeature_RNA', 'percent.mt',
						'age_group_adult','age_group_early_dev','age_group_late_dev')
	} else {
		features <- c(grep('scVI', colnames(embeddings), value = TRUE),
						'UMAP_1', 'UMAP_2', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')
	}
	write(features, write_features_file)

	# train on pre-labelled cells
	pickle <- paste0(temp_name, '_', rand_num, '.pickle' )
	system(glue('{conda_dir}/envs/integrate_scRNA/bin/python3.6 {git_dir}/src/cell_type_predictor.py train --predProbThresh {probThresh} --workingDir {working_dir} --inputMatrix  {embeddings_file}  --trainedModelFile  {pickle}  --featureCols  {write_features_file}'))

	# apply model to predict labels for all cells
	predictions_file <- temp_name
	system(glue('{conda_dir}/envs/integrate_scRNA/bin/python3.6 {git_dir}/src/cell_type_predictor.py predict --predProbThresh  {probThresh} --workingDir {working_dir} --inputMatrix  {full_embeddings_file} --trainedModelFile  {pickle} --predictions  {predictions_file}'))

	# import in predictions
	predictions <- read_tsv(predictions_file)

	predictions
}

# full data set
predictions <- run_xgboost_py(out, full_out, 'fullTemp', tm = FALSE, probThresh = 0.5)
umapX <- left_join(umap, predictions %>% select(Barcode, CellType_predict = CellType), by = 'Barcode')
## remove tabula muris
umapX <- umapX %>% filter(study_accession != 'SRP131661')

# just tabula muris to fill in missing "TabulaMurisCellType"
predictions_tm <- run_xgboost_py(out_tm, full_tm, 'tmTemp', tm = TRUE, probThresh = 0.5)
umapX2 <- left_join(umap %>% filter(study_accession == 'SRP131661'), predictions_tm %>% select(Barcode, TabulaMurisCellType_predict = CellType), by = 'Barcode')

# glue together
umapX3 <- bind_rows(umapX, umapX2)
if (nrow(umapX3) != nrow(umap)){
	print('Cells lost!!!')
	stop()
}

umap <- umapX3 
umap$CellType_predict <- gsub('None', NA, umap$CellType_predict)
save(predictions, predictions_tm, file = args[3])
save(umap, file = args[4])
