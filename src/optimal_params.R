library(tidyverse)
library(jsonlite)
#args <- '/data/swamyvs/scEiaD/rson_tmp/ji_whfuw.json'
args <- commandArgs(trailingOnly = TRUE)
rule <- read_json(args[1])
partitions_run = str_split(rule$input$perf_data %>% unlist, '__') %>% 
  unlist %>% 
  .[grepl('partition', .)] %>% 
  str_split('-') %>% sapply(function(x) x[2]) %>% unique()
#library(cowplot)
#library(splancs)
# load annotations
# load('cell_info_labelled.Rdata')
# calculate of matrix of points
area_chull <- function(x,y){
  matrix <- cbind(x,y) %>% as.matrix()
  ch = chull(matrix)
  coords <- matrix[c(ch, ch[1]),]
  Polygon(coords)@area
}

perf_all <- list()
for (x in partitions_run){
  #for (x in c('.*TabulaDroplet.*', '.*onlyWELL.*', '.*universe.*')){
  print(x)
  # load all clustering params -----
  files <- list.files('pipeline_data/perf_metrics',
                      pattern = paste0(x, '.*Rdata'),
                      full.names = TRUE)
  for (i in files){
    print(i)
    load(i)
    if (x == 'onlyWELL') {
      table_score <- rbind(c('LISI', scores$LISI_batch %>% pull(1) %>% mean(), 'Batch'),
                           c('LISI', scores$LISI_cluster %>% pull(1) %>% mean(), 'Cluster'),
                           c('Silhouette', scores$silhouette_batch, 'Batch'),
                           c('Silhouette', scores$silhouette_cluster, 'Cluster'),
                           scores$RI %>% enframe(name = 'score') %>% mutate('Cell Type' = 'CellType-Cluster'))
    } else {
      table_score <- rbind(c('LISI', scores$LISI_batch %>% pull(1) %>% mean(), 'Batch'),
                           c('LISI', scores$LISI_celltype %>% pull(1) %>% mean(), 'CellType'),
                           c('LISI', scores$LISI_cluster %>% pull(1) %>% mean(), 'Cluster'),
                           c('LISI', scores$LISI_subcelltype %>% pull(1) %>% mean(), 'SubCellType'),
                           c('Silhouette', scores$silhouette_batch, 'Batch'),
                           c('Silhouette', scores$silhouette_celltype, 'CellType'),
                           c('Silhouette', scores$silhouette_cluster, 'Cluster'),
                           c('Silhouette', scores$silhouette_subcelltype, 'SubCellType'),
                           scores$RI %>% enframe(name = 'score') %>% mutate('Cell Type' = 'CellType-Cluster'))
    }
    colnames(table_score) <- c('Score', 'Value', 'Group')
    table_score$Value <- as.numeric(as.character(table_score$Value))
    
    norm = str_extract(i, 'transform-[a-zA-Z]+') %>% gsub('transform-','',.)
    nf = str_extract(i, 'n_features-\\d+') %>% gsub('n_features-', '', .) %>% as.numeric()
    dims = str_extract(i, 'dims-\\d+') %>% gsub('dims-', '', .) %>% as.numeric()
    method = str_extract(i, 'method-[a-zA-Z]+') %>% gsub('method-','',.)
    knn = str_extract(i, 'knn-\\d+\\.\\d+|knn-\\d+') %>% gsub('knn-', '', .) %>% as.numeric()
    epochs = str_extract(i, 'epochs-\\d+') %>% gsub('epochs-','',.) %>% as.numeric()
    table_score$clusterN = max(as.numeric(scores$cluster_count$Cluster))
	table_score$clusterMedian = median(as.numeric(scores$cluster_count$Count))
	table_score$dims = dims
    table_score$nf = nf
    table_score$knn <- knn
    table_score$epochs <- epochs
    table_score$method <- method	
    table_score$normalization <- norm
    table_score$subset <- x
    table_score$set <- x
    perf_all[[i]] <- table_score
  }
  #perf_all <- perf_all %>% bind_rows()
}
perf_one <- perf_all %>% bind_rows() %>% unique()


# scIB
files = rule$input$scib_data
data <- list()
for (i in files){
  norm = str_extract(i, 'transform-[a-zA-Z]+') %>% gsub('transform-','',.)
  nf = str_extract(i, 'n_features-\\d+') %>% gsub('n_features-', '', .) %>% as.numeric()
  dims = str_extract(i, 'dims-\\d+') %>% gsub('dims-', '', .) %>% as.numeric()
  method = str_extract(i, 'method-[a-zA-Z]+') %>% gsub('method-','',.)
  knn = str_extract(i, 'knn-\\d+\\.\\d+|knn-\\d+') %>% gsub('knn-', '', .) %>% as.numeric() 
  epochs = str_extract(i, 'epochs-\\d+') %>% gsub('epochs-','',.) %>% as.numeric()
  set = str_split(i %>% unlist, '__') %>% 
    unlist %>% 
    .[grepl('partition', .)] %>% 
    str_split('-') %>% sapply(function(x) x[2])
  
  table_score <- read_csv(i, col_names = FALSE)
  table_score$dims = dims
  table_score$nf = nf
  table_score$knn <- knn
  table_score$method <- method
  table_score$normalization <- norm
  table_score$epochs <- epochs
  table_score$set <- set
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

perf <- bind_rows(perf_one %>% filter(Score != 'ARI'), perf_two) %>% arrange(Score, desc(Value))
save(perf, file = rule$output$merged_stats)
