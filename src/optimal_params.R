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
x2scf = '.*n_features2000.*scran.*fast'
x2stf = '.*n_features2000.*standard.*fast'
x2cs = '.*n_features2000.*count.*scVI'
x5 = '.*n_features5000.*count.*scVI'
x10 = '.*n_features10000.*count.*scVI'

umap_mega <- list()
for (x in c(x1, x2cs, x5, x10)){
  print(x)
  # load all clustering params -----
  files <- list.files('cluster/', 
                      pattern = paste0(x, '.*.cluster.Rdata'), 
                      full.names = TRUE)
  cluster_all <- list()
  for (i in files){
    load(i)
    colnames(meta) <- c('Barcode', 'cluster')
    cluster_all[[i]] <- meta
    nf = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
    knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
    cluster_all[[i]]$dims = dims
    cluster_all[[i]]$knn = knn
    cluster_all[[i]]$nf = nf
    cluster_all[[i]]$cluster_sum = max(meta$cluster)
    cluster_all[[i]]$median_cluster_n = max(meta$cluster)
    cluster_all[[i]]$method <- method
  }
  cluster_all <- cluster_all %>% bind_rows()
  
  
  
  # load umap file with one cluster param -----
  files <- list.files('umap/', 
                      pattern =  paste0(x, '.*mindist0.3.*'), 
                      full.names = TRUE)
  umap_all <- list()
  count = 1
  for (i in files){
    print(count)
    count = count + 1
    load(i)
    nn = str_extract(i, 'nneighbors\\d+') %>% gsub('nneighbors', '', .) %>% as.numeric()
    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
    dist = str_extract(i, 'mindist0.\\d+') %>% gsub('mindist', '', .) %>% as.numeric()
    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
    norm = str_extract(i, 'n_features\\d+__[^\\W_]+') %>% gsub('n_features\\d+__','',.)
    umap_all[[i]] <- umap
    umap_all[[i]]$mindist = dist
    umap_all[[i]]$nneighbors = nn
    umap_all[[i]]$dims = dims
    umap_all[[i]]$method = method
    umap_all[[i]]$normalization = norm
    umap_all[[i]]$total_area = area_chull(umap_all[[i]]$UMAP_1, umap_all[[i]]$UMAP_2)
    d <- dims
    for (k in c(cluster_all$knn %>% unique())){
      umap_all[[paste0(i,'_', k)]] <- umap_all[[i]] %>% 
        select(-cluster) %>% 
        left_join(., cluster_all %>% filter(dims == d, knn == k), by = 'Barcode') %>% 
        filter(!is.na(CellType), !is.na(cluster)) %>%
        mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType)) %>%
        group_by(total_area, cluster, CellType, organism, study_accession) %>% 
        summarise(Count = n(),
                  x = mean(UMAP_1), 
                  y = mean(UMAP_2),
                  area = area_chull(UMAP_1, UMAP_2),
                  mindist = unique(mindist),
                  nneighbors = unique(nneighbors),
                  dims = unique(dims)) %>% 
        filter(Count > 10) %>% 
        summarise(study_count = n(), 
                  Count = sum(Count),
                  x = mean(x), 
                  y = mean(y),
                  area = mean(area),
                  mindist = unique(mindist),
                  nneighbors = unique(nneighbors),
                  dims = unique(dims)) %>% 
        summarise(organism_count = n(), 
                  Count = sum(Count), 
                  study_count = sum(study_count),
                  x = mean(x), 
                  y = mean(y),
                  area = mean(area),
                  mindist = unique(mindist),
                  nneighbors = unique(nneighbors),
                  dims = unique(dims)) %>% 
        mutate(freq = Count / sum(Count),
               area_scaled = area / total_area) %>% 
        ungroup() %>% 
        group_by(dims, cluster) %>% 
        top_n(n = 1, wt = freq) 
      umap_all[[paste0(i,'_', k)]]$knn <- k
      umap_all[[paste0(i,'_', k)]]$method <- method
      umap_all[[paste0(i,'_', k)]]$normalization <- norm
    }
    umap_all[[i]] <- NULL
  }
  umap_one <- umap_all %>% bind_rows()
  umap_one$nfeatures <- x
  umap_mega[[paste0(x, '_', method, '_', norm)]] <- umap_one
  
}
umap_one <- umap_mega %>% bind_rows()# %>% filter(dims != 8)
save(umap_one, file = 'umap_one.Rdata')


# process ari / silhouette / etc data 
perf_all <- list()
for (x in c(x1, x2cs, x5, x10)){
  print(x)
  # load all clustering params -----
  files <- list.files('perf_metrics',
                      pattern = paste0(x, '.*Rdata'),
                      full.names = TRUE)
  for (i in files){
    load(i)
	names(scores$silhouette_CellType) <-  scores$umap_cutdown %>% group_by(CellType) %>% count() %>% pull(1)
    names(scores$RI_species) <-  scores$umap_cutdown %>% group_by(organism) %>% count() %>% pull(1)
	lisi <- cbind(scores$LISI, scores$umap_cutdown$CellType) %>% group_by(`scores$umap_cutdown$CellType`) %>% summarise(value = mean(CellType)) %>% mutate(score = 'LISI')
	colnames(lisi)[1] <- 'Cell Type'	
	lisi[,1] <- lisi[,1] %>% pull(1) %>% as.character()
	table_score <- rbind(lisi, c('All', mean(lisi$value), 'LISI')) %>% mutate(value = as.numeric(value))
    table_score <- bind_rows(table_score, 
							 scores$RI %>% enframe(name = 'score') %>% mutate('Cell Type' = 'All'))
	silhouette <- scores$silhouette_CellType %>% bind_rows() %>% t() %>% data.frame()
	silhouette$`Cell Type` = row.names(silhouette)
	silhouette$score <- 'Silhouette'
	colnames(silhouette)[1] <- 'value'
	silhouette <- as_tibble(silhouette)
	table_score <- bind_rows(table_score, silhouette)
	table_score <- rbind(table_score, c('All', scores$silhouette, 'Silhouette'))

	nf = str_extract(i, 'n_features\\d+') %>% gsub('n_features', '', .) %>% as.numeric()
    dims = str_extract(i, 'dims\\d+') %>% gsub('dims', '', .) %>% as.numeric()
    #knn = str_extract(i, 'knn\\d+') %>% gsub('knn', '', .) %>% as.numeric()
    method = str_extract(i, 'batch__[^\\W_]+') %>% gsub('batch__','',.)
    table_score$dims = dims
    #table_score$knn = knn
    table_score$nf = nf
	table_score$method <- method	

	perf_all[[i]] <- table_score
  }
  #perf_all <- perf_all %>% bind_rows()
}
