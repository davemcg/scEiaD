library(tidyverse)
library(cowplot)
library(splancs)
# load annotations
# load('cell_info_labelled.Rdata')
# calculate of matrix of points
area_chull <- function(x,y){
  matrix <- cbind(x,y) %>% as.matrix()
  ch = chull(matrix)
  coords <- matrix[c(ch, ch[1]),]
  Polygon(coords)@area
}
# regex to pick projection types
x1 <- '.*n_features1000.*count.*scVI'
x2scf = '.*n_features2000.*scran'
x2stf = '.*n_features2000.*standard'
x2cs = '.*n_features2000.*count.*scVI'
x5 = '.*n_features5000.*count.*scVI'
x10 = '.*n_features10000.*count.*scVI'

#umap_mega <- list()
#for (x in c('.*onlyDROPLET.*')){
#  print('cluster data')
#  # load all clustering params -----
#  files <- list.files('cluster/', 
#                      pattern = paste0(x, '.*.cluster.Rdata'), 
#                      full.names = TRUE)
#  cluster_all <- list()
#  count = 1
#  for (i in files){
#    print(count)
#    count = count + 1
#    load(i)
#    colnames(meta) <- c('Barcode', 'cluster')
#    cluster_all[[i]] <- meta
#    nf = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
#    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
#    knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
#    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
#    norm = str_extract(i, 'n_features\\d+__[^\\W_]+') %>% gsub('n_features\\d+__','',.)
#	meta$cluster <- meta$cluster %>% as.character() %>% as.numeric()
#	cluster_all[[i]]$dims = dims
#    cluster_all[[i]]$knn = knn
#    cluster_all[[i]]$nf = nf
#    cluster_all[[i]]$cluster <- meta$cluster
#	cluster_all[[i]]$cluster_sum = max(meta$cluster)
#    cluster_all[[i]]$median_cluster_n = max(meta$cluster)
#    cluster_all[[i]]$method <- method
#	cluster_all[[i]]$norm <- norm
#  }
#  
#  cluster_all <- cluster_all %>% bind_rows() %>% as_tibble(.name_repair = 'unique')
#  cluster_all <- cluster_all %>% select(-contains('...')) 
#  # load umap file with one cluster param -----
#  files <- list.files('umap/', 
#                      pattern =  paste0(x, '.*mindist0.3.*'), 
#                      full.names = TRUE)
#  umap_all <- list()
#  count = 1
#  for (i in files){
#    print('umap time')
#    print(count)
#    count = count + 1
#    load(i)
#    if (!"cluster" %in% colnames(umap)){
#    	umap$cluster <- umap$clusters
#	}
#    nn = str_extract(i, 'nneighbors\\d+') %>% gsub('nneighbors', '', .) %>% as.numeric()
#    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
#    dist = str_extract(i, 'mindist0.\\d+') %>% gsub('mindist', '', .) %>% as.numeric()
#    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
#    norm = str_extract(i, 'n_features\\d+__[^\\W_]+') %>% gsub('n_features\\d+__','',.)
#	nfeatures = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
#    umap_all[[i]] <- umap
#    umap_all[[i]]$mindist = dist
#    umap_all[[i]]$nneighbors = nn
#    umap_all[[i]]$dims = dims
#    umap_all[[i]]$method = method
#    umap_all[[i]]$normalization = norm
#    umap_all[[i]]$total_area = area_chull(umap_all[[i]]$UMAP_1, umap_all[[i]]$UMAP_2)
#    umap_all[[i]]$nfeatures <- nfeatures
#    d <- dims
#	m <- method
#    n <- norm
#    for (k in c(cluster_all %>% filter(dims == d,  method == m, nf == nfeatures, norm == n) %>% pull(knn) %>% unique())){
#      umap_all[[paste0(i,'_', k)]] <- umap_all[[i]] %>% 
#        select(-cluster) %>% 
#		left_join(., cluster_all %>% filter(dims == d, knn == k, method == m, nf == nfeatures, norm == n), by = 'Barcode') %>% 
#        filter(!is.na(CellType), !is.na(cluster)) %>%
#        filter(!CellType %in% c('Astrocytes', 'Red Blood Cells', 'Doublet', 'Doublets',  'Fibroblasts')) %>% 
#		mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>%
#        group_by(total_area, cluster, CellType, organism, study_accession) %>% 
#        summarise(Count = n(),
#                  x = mean(UMAP_1), 
#                  y = mean(UMAP_2),
#                  area = area_chull(UMAP_1, UMAP_2),
#                  mindist = unique(mindist),
#				  normalization = unique(normalization), 
#                  method = unique(method),
#				  nneighbors = unique(nneighbors),
#				  nfeatures = unique(nfeatures),
#                  dims = unique(dims)) %>% 
#        filter(Count > 10) %>% 
#        summarise(study_count = n(), 
#                  Count = sum(Count),
#                  x = mean(x), 
#                  y = mean(y),
#                  area = mean(area),
#                  mindist = unique(mindist),
#				  normalization = unique(normalization),
#                  method = unique(method),
#                  nneighbors = unique(nneighbors),
#				  nfeatures = unique(nfeatures),
#                  dims = unique(dims)) %>% 
#        summarise(organism_count = n(), 
#                  Count = sum(Count), 
#                  study_count = sum(study_count),
#                  x = mean(x), 
#                  y = mean(y),
#                  area = mean(area),
#                  mindist = unique(mindist),
#				  normalization = unique(normalization),
#                  method = unique(method),
#                  nneighbors = unique(nneighbors),
# 				  nfeatures = unique(nfeatures),
#                  dims = unique(dims)) %>% 
#        mutate(freq = Count / sum(Count),
#               area_scaled = area / total_area) %>% 
#        ungroup() %>% 
#        group_by(dims, cluster) %>% 
#        top_n(n = 1, wt = freq) 
#      umap_all[[paste0(i,'_', k)]]$knn <- k
#      #umap_all[[paste0(i,'_', k)]]$method <- method
#      #umap_all[[paste0(i,'_', k)]]$normalization <- norm
#    }
#    umap_all[[i]] <- NULL
#  }
#  umap_one <- umap_all %>% bind_rows()
#  umap_mega[[paste0(x, '_', method, '_', norm)]] <- umap_one
#}
#umap_one <- umap_mega %>% bind_rows()# %>% filter(dims != 8)
#save(umap_one, file = 'umap_one.Rdata')


# process ari / silhouette / etc data 
perf_all <- list()
for (x in c('.*onlyDROPLET.*')){
  print(x)
  # load all clustering params -----
  files <- list.files('perf_metrics',
                      pattern = paste0(x, '.*Rdata'),
                      full.names = TRUE)
  for (i in files){
    load(i)
	table_score <- rbind(c('LISI', scores$LISI_batch %>% pull(1) %>% mean(), 'Batch'),
								c('LISI', scores$LISI_celltype %>% pull(1) %>% mean(), 'CellType'),
								c('LISI', scores$LISI_cluster %>% pull(1) %>% mean(), 'Cluster'),
								c('Silhouette', scores$silhouette_batch, 'Batch'),
								c('Silhouette', scores$silhouette_celltype, 'CellType'),
								c('Silhouette', scores$silhouette_cluster, 'Cluster'),
								scores$RI %>% enframe(name = 'score') %>% mutate('Cell Type' = 'CellType-Cluster'))
							
	colnames(table_score) <- c('Score', 'Value', 'Group')
	table_score$Value <- as.numeric(as.character(table_score$Value))


	norm = str_extract(i, 'n_features\\d+__[^\\W_]+') %>% gsub('n_features\\d+__','',.)
	nf = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
    knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
	table_score$dims = dims
    table_score$nf = nf
	table_score$knn <- knn
	table_score$method <- method	
	table_score$normalization <- norm

	perf_all[[i]] <- table_score
  }
  #perf_all <- perf_all %>% bind_rows()
}
perf_one <- perf_all %>% bind_rows() %>% unique()


# scIB
files = list.files('scIB_stats/', pattern = 'csv', full.names = TRUE)
data <- list()
for (i in files){
  	norm = str_extract(i, 'n_features\\d+__[^\\W_]+') %>% gsub('n_features\\d+__','',.)
    nf = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
    knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
	table_score <- read_csv(i, col_names = FALSE)
	table_score$dims = dims
    table_score$nf = nf
    table_score$knn <- knn
    table_score$method <- method
    table_score$normalization <- norm
	colnames(table_score)[1:2] <- c('Score', 'Value')
	data[[i]] <- table_score
}
perf_two <- data %>% bind_rows() %>%  
			mutate(Group = case_when(Score == 'nmi' ~ 'CellType-Cluster',
									 Score == 'nmi_sub' ~ 'SubCellType-Cluster',
									 Score == 'ari' ~ 'CellType-Cluster',
							         Score == 'ari_sub' ~ 'SubCellType-Cluster',
									Score == 'pcr' ~ 'After-Before')) %>%
			mutate(Score = gsub('_sub','',Score)) %>% 
			mutate(Score = toupper(Score))

perf <- bind_rows(perf_one %>% filter(Score != 'ARI'), perf_two)
save(perf, file = 'metrics_onlyDROPLET_2020_06_02.Rdata')
