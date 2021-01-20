library(tidyverse)
library(glue)
library(reticulate)
args <- commandArgs(trailingOnly = TRUE)

load(args[1]) # seurat obj post integration
load(args[2]) # umap metadata
label_id_col = args[5]# label_id_col is the name of a column that will be created and used interally within the python script,
label_name_col = args[6]# label_name_col must a be a column in embedding/umap file - CellType, cluster etc 
model_outfile <- args[7]


conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
working_dir = Sys.getenv('SCIAD_WORKING_DIR')

out <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode') %>% left_join(umap, by = 'Barcode')
out_tm <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode') %>% left_join(umap, by = 'Barcode') %>% filter(study_accession == 'SRP131661')
out_tm$CellType <- out_tm$TabulaMurisCellType
remove_low <- out_tm %>% group_by(CellType) %>% dplyr::count() %>% filter(n < 40) %>% pull(CellType)
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
		cutoff = sd(D) * 3 
		bc_retain[[i]] <- D[D < cutoff] %>% names()
	}
	out %>% filter(Barcode %in% (bc_retain %>% unlist()))
}
#out$id <- out$CellType %>% as.factor() %>% as.numeric()
#out$id <- out_tm$CellType %>% as.factor() %>% as.numeric()

full_out <- out
full_tm <- out_tm
out <- rm_outlier(out)
out_tm <- rm_outlier(out_tm)

run_xgboost_py <- function(embeddings, full_embeddings, model_outfile, tm = FALSE, probThresh = 0.8, 
                           label_id_col=label_id_col, label_name_col=label_name_col){
  
	if (!tm) {
		features <- c(grep('scVI', colnames(embeddings), value = TRUE),
						'UMAP_1', 'UMAP_2', 'nCount_RNA', 'nFeature_RNA', 'percent.mt',
						'age_group_adult','age_group_early_dev','age_group_late_dev')
	} else {
		features <- c(grep('scVI', colnames(embeddings), value = TRUE),
						'UMAP_1', 'UMAP_2', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')
	}
	#write(features, write_features_file)
	
	# train on pre-labelled cells
	pickle <- paste0(model_outfile, '_.pickle' )
  use_python(glue('{conda_dir}/envs/integrate_scRNA/bin/python'))
  source_python( glue('{git_dir}/src/cell_type_predictor.py'))
  train_test_predictions <-  scEiaD_classifier_train(inputMatrix=embeddings, labelIdCol=label_id_col, labelNameCol=label_name_col,  trainedModelFile=pickle,
        featureCols=features, predProbThresh=probThresh, generateProb=TRUE)
  
	# apply model to predict labels for all cells
	
	full_embedding_predictions <-scEiaD_classifier_predict(inputMatrix=full_embeddings,labelIdCol=label_id_col, labelNameCol=label_name_col, trainedModelFile=pickle,
	                                                       featureCols=features,  predProbThresh=probThresh)
	out <- list()
	out[['predictions']] <- full_embedding_predictions
	out[['test_data']] <- train_test_predictions$test_probs_df
	out[['training_data']] <- train_test_predictions$train_probs_df
	out
}

# full data set

model_out <- run_xgboost_py(out, full_out, 'fullTemp', tm = FALSE, probThresh = 0.5)
predictions <- model_out$predictions
umapX <- left_join(umap, predictions %>% select(Barcode, CellType_predict = CellType), by = 'Barcode')
## remove tabula muris
umapX <- umapX %>% filter(study_accession != 'SRP131661')

# just tabula muris to fill in missing "TabulaMurisCellType"
modelTM_out <- run_xgboost_py(out_tm, full_tm, 'tmTemp', tm = TRUE, probThresh = 0.5)
predictions_tm <- modelTM_out$predictions
umapX2 <- left_join(umap %>% filter(study_accession == 'SRP131661'), predictions_tm %>% select(Barcode, TabulaMurisCellType_predict = CellType), by = 'Barcode')

# glue together
umapX3 <- bind_rows(umapX, umapX2)
if (nrow(umapX3) != nrow(umap)){
	print('Cells lost!!!')
	stop()
}

umap <- umapX3 
umap$CellType_predict <- gsub('None', NA, umap$CellType_predict)
save(predictions, predictions_tm, model_out, modelTM_out, file = args[3])
save(umap, file = args[4])
