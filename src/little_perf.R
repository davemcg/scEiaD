library(aricode)
library(Seurat)
library(tidyverse)
library(cluster)

args = commandArgs(trailingOnly=TRUE)

load(args[1])


sample <- str_extract(args[1], '\\/.*__nf') %>% gsub('/|__nf','',.)
dims <-  str_extract(args[1], 'dims\\d+') %>% gsub('dims','',.)
norm <- str_extract(args[1], 'norm\\w+') %>% gsub('norm','',.)

meta <- seurat@meta.data %>% as_tibble()

ari_CellType <- try({ ARI(meta %>% filter(!is.na(CellType)) %>% pull(CellType), meta %>% filter(!is.na(CellType)) %>% pull(seurat_clusters)) })
ari_CellType_predict <- try( { ARI(meta %>% filter(!is.na(CellType_predict)) %>% pull(CellType_predict), meta %>% filter(!is.na(CellType_predict)) %>% pull(seurat_clusters)) })

nmi_CellType <- try({ NMI(meta %>% filter(!is.na(CellType)) %>% pull(CellType), meta %>% filter(!is.na(CellType)) %>% pull(seurat_clusters)) })
nmi_CellType_predict <- try({ NMI(meta %>% filter(!is.na(CellType_predict)) %>% pull(CellType_predict), meta %>% filter(!is.na(CellType_predict)) %>% pull(seurat_clusters)) })

# silhouette
asw <- function(seurat_obj, meta_column, dims){

	bc_rm <- seurat_obj@meta.data[,meta_column]
	bc_rm <- is.na(bc_rm)
	seuratCT <- seurat_obj[, !bc_rm]
	dist.matrix <- dist(x = Embeddings(object = seuratCT[['pca']])[, 1:dims])
	clusters <- seuratCT@meta.data[,meta_column]
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	sil[,3] %>% mean()
}
asw_CT <- try({asw(seurat, 'CellType', dims)})
asw_CTpred <- try({asw(seurat, 'CellType_predict', dims)})
if (class(asw_CT) == 'try-error'){
 asw_CT = 0
}
if (class(asw_CTpred) == 'try-error'){
 asw_CTpred = 0
}

table <- as_tibble(c(ari_CellType, ari_CellType_predict, nmi_CellType, nmi_CellType_predict, asw_CT, asw_CTpred))
table$test <- c('ARI','ARI','NMI','NMI', 'ASW','ASW')
table$comparison <- c('CellType - Cluster','CellType_predict - Cluster','CellType - Cluster','CellType_predict - Cluster', 'CellType','CellType_predict')
table$sample <- sample
table$dims <- dims
table$norm <- norm

write_tsv(table, file = args[2])

Sys.sleep(300)
