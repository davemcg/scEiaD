library(reticulate)
library(matrixStats)
args <- commandArgs(trailingOnly = TRUE)
library(glue)
conda_dir = Sys.getenv('SCIAD_CONDA_DIR')
git_dir = Sys.getenv('SCIAD_GIT_DIR')
working_dir = Sys.getenv('SCIAD_WORKING_DIR')
use_python(glue('{conda_dir}/envs/integrate_scRNA/bin/python'), required = TRUE)
source_python( glue('{git_dir}/src/cell_type_predictor.py'))
library(tidyverse)
library(lsa) # cosine

load(args[1]) # seurat obj post integration
load(args[2]) # umap metadata
label__id__col = args[5]# label_id_col is the name of a column that will be created and used interally within the python script,
label__name__col = args[6]# label_name_col must a be a column in embedding/umap file - CellType, cluster etc 
model_outfile <- args[7]


print(args)

out <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode') %>% left_join(umap %>% mutate(CellType = gsub("Cone Bipolar Cells", "Bipolar Cells", CellType)), by = 'Barcode')
out_tm <- integrated_obj@reductions$scVI@cell.embeddings %>% as_tibble(rownames = 'Barcode') %>% left_join(umap %>% mutate(CellType = gsub("Cone Bipolar Cells", "Bipolar Cells", CellType)), by = 'Barcode') %>% filter(study_accession == 'SRP131661')
out_tm$CellType <- out_tm$TabulaMurisCellType

# ensure Age is numeric
out$Age <- as.numeric(out$Age)
out_tm$Age <- as.numeric(out_tm$Age)

# cluster outlier
# removes cell tyoe calls that are very rare in a cluster
cluster_outlier <- function(out){
	groupings <- out %>% 
					filter(!is.na(CellType)) %>%  
					group_by(cluster, CellType) %>% 
					summarise(Count = n()) %>% 
					mutate(Perc = Count / sum(Count)) %>% 
					filter(Perc > 0.025)
	out %>% 
			left_join(groupings %>% 
			select(cluster, CellTypeKEEP = CellType)) %>% 
			filter(CellType == CellTypeKEEP) %>% 
			select(-CellTypeKEEP)
}

# remove celltype outliers
# we will removve 5% furthest from cosine dist center
rm_outlier <- function(out, TM = FALSE){
	bc_retain <- list()
	for (i in out$CellType %>% unique){
		print(i)
		temp <- out %>% filter(CellType == i) %>% select(contains('scVI_')) %>% as.matrix()
		row.names(temp) <- out %>% filter(CellType == i) %>% pull(Barcode)
		scVI_mean <- colMeans(temp)
		#euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
		D <- apply(as.matrix(temp), 1, function(x) cosine(scVI_mean, x)) #cosine dist
		cutoff = quantile(D, probs = seq(0,1,0.01))[96]
		#print(glue("Retaining {D[D < cutoff] %>% names() %>% length()}" )) 
		#print(glue("Removing {D[D > cutoff] %>% names() %>% length()}" ))
		bc_retain[[i]] <- D[D < cutoff] %>% names()
	}
	out %>% filter(Barcode %in% (bc_retain %>% unlist()))
}
#out$id <- out$CellType %>% as.factor() %>% as.numeric()
#out$id <- out_tm$CellType %>% as.factor() %>% as.numeric()

# remove rare celltypes

full_out <- out
full_tm <- out_tm
out <- cluster_outlier(out) %>% rm_outlier(.)
#out <- rm_outlier(out)
out_tm <- rm_outlier(out_tm)

# remove low n celltypes
remove_low <- out %>% group_by(CellType) %>% dplyr::count() %>% filter(n < 20) %>% pull(CellType)
out <- out %>% filter(!CellType %in% remove_low)

remove_low <- out_tm %>% group_by(CellType) %>% dplyr::count() %>% filter(n < 20) %>% pull(CellType)
out_tm <- out_tm %>% filter(!CellType %in% remove_low)

run_xgboost_py <- function(embeddings, full_embeddings, model_outfile, tm = FALSE, probThresh = 0.8, 
                           label_id_col=label__id__col, label_name_col=label__name__col){
  
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
	pickle <- #paste0(model_outfile, '.pickle' )
	pickle <- model_outfile
  train_test_predictions <-  scEiaD_classifier_train(inputMatrix=embeddings, labelIdCol=label_id_col, labelNameCol=label_name_col,  trainedModelFile=pickle,
        featureCols=features, predProbThresh=probThresh, generateProb=TRUE,
         bad_cell_types = list("RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cell", "Droplet", "Droplets", 'Doublet', 'Doublets', 'Choriocapillaris','Artery'))
  
	# apply model to predict labels for all cells
	
	full_embedding_predictions <-scEiaD_classifier_predict(inputMatrix=full_embeddings,labelIdCol=label_id_col, labelNameCol=label_name_col, trainedModelFile=pickle,
	                                                       featureCols=features,  predProbThresh=probThresh)
	out <- list()
	out[['predictions']] <- full_embedding_predictions
	out[['test_data']] <- train_test_predictions$test_probs_df %>% py_to_r()
	out[['training_data']] <- train_test_predictions$train_probs_df %>% py_to_r()
	out
}

# full data set

model_out <- run_xgboost_py(out, full_out, model_outfile, tm = FALSE, probThresh = 0.5)
predictions <- model_out$predictions
predictions$CellType_predict_max_prob <- model_out$predictions %>% select_if(is.numeric) %>% select(-CellTypeID) %>% mutate(CellType_predict_prob = do.call(pmax, select_if(., is.numeric))) %>% pull(CellType_predict_prob)
umapX <- left_join(umap, predictions %>% select(Barcode, CellType_predict = CellType, CellType_predict_max_prob), by = 'Barcode')
## remove tabula muris
umapX <- umapX %>% filter(study_accession != 'SRP131661')

# just tabula muris to fill in missing "TabulaMurisCellType"
modelTM_out <- run_xgboost_py(out_tm, full_tm, paste0(model_outfile, 'TEMP'), tm = TRUE, probThresh = 0.5)
predictions_tm <- modelTM_out$predictions
predictions_tm$CellType_predict_max_prob <- modelTM_out$predictions %>% select_if(is.numeric) %>% select(-CellTypeID) %>% mutate(CellType_predict_prob = do.call(pmax, select_if(., is.numeric))) %>% pull(CellType_predict_prob)
umapX2 <- left_join(umap %>% filter(study_accession == 'SRP131661'), predictions_tm %>% select(Barcode, TabulaMurisCellType_predict = CellType, TabulaMurisCellType_predict_max_prob = CellType_predict_max_prob), by = 'Barcode')

# glue together
umapX3 <- bind_rows(umapX, umapX2)
if (nrow(umapX3) != nrow(umap)){
	print('Cells lost!!!')
	stop()
}

umap <- umapX3 
umap$CellType_predict <- gsub('None', NA, umap$CellType_predict)

# quick accuracy
cell_type2id <- predictions %>% select(CellTypeID, CellType) %>% distinct %>% filter(CellType!='None')
target_cell_types <-unique(cell_type2id$CellType)
bc_col <- grep('Barcode', colnames(predictions))
test_predictions <- predictions %>%
  mutate(max_pred_prob = rowMaxs(.[,-(bc_col:ncol(predictions))] %>% as.matrix )) %>%
  select(-CellTypeID) %>%
  rename(PredCellType = CellType) %>%
  inner_join(umap %>% 
		mutate(CellType = gsub("Cone Bipolar Cells", "Bipolar Cells", CellType)) %>%
		select(Barcode, TrueCellType = CellType)) %>%
  filter(!is.na(TrueCellType)) %>%
  mutate(pred_correct  = ifelse(PredCellType == TrueCellType, 'Correct', 'incorrect')) %>% 
  filter(PredCellType != 'None')

accuracy <- test_predictions %>% ungroup() %>% group_by(TrueCellType) %>% summarise(score = sum(pred_correct == 'Correct') / n()) %>% select(TrueCellType, score)

save(predictions, predictions_tm, model_out, modelTM_out, accuracy, file = args[3])
save(umap, file = args[4])
