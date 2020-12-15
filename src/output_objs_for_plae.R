library(pool)

library(tidyverse)
library(Seurat)
args = commandArgs(trailingOnly=TRUE)

# get umap embeddings
load(args[1])
#load('seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umap.Rdata')
temp <- integrated_obj
load(args[2])
#load('seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter.seuratV3.Rdata')

#integrated_obj@reductions$scVI <- CreateDimReducObject(embeddings = temp@reductions$scVI@cell.embeddings,
#                                                  loadings = temp@reductions$scVI@feature.loadings,
#                                                  projected = temp@reductions$scVI@feature.loadings.projected,
#                                                  key = 'scVI_',
#                                                  assay = 'RNA')
integrated_obj@reductions$scviUMAP <- CreateDimReducObject(embeddings = temp@reductions$scviUMAP@cell.embeddings,
                                                      loadings = temp@reductions$scviUMAP@feature.loadings,
                                                      projected = temp@reductions$scviUMAP@feature.loadings.projected,
                                                      key = 'scVIUMAP_',
                                                      assay = 'RNA')
load(args[3])
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umap.Rdata')

new_cols <- colnames(umap)[!(colnames(umap) %in% colnames(integrated_obj@meta.data))]

for (i in new_cols){ 
	new_data <- umap[,i] %>% pull(1)
	names(new_data) <- umap$Barcode
	integrated_obj <- AddMetaData(integrated_obj, new_data, col.name = i)
}

load(args[4])
#load('cluster/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__knn7.cluster.Rdata')

integrated_obj <- AddMetaData(integrated_obj, meta[,2] %>% pull(1), col.name = 'cluster')

scEiaD_droplet <- integrated_obj

save(scEiaD_droplet, file = 'site/scEiaD_all_seurat_v3.Rdata')

#load(args[5])
#load('well_data_seurat_obj_labelled.Rdata')
#seurat_obj <- AddMetaData(seurat_obj, seurat_obj@meta.data$seurat_clusters, col.name=  'cluster')
#scEiaD_well <- seurat_obj
#save(scEiaD_well, file = 'site/scEiaD_well_seurat_v3.Rdata')

# rename raw data
library(Matrix)
load(args[5])
#load('seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features1000__counts__raw__batch__preFilter.seuratV3.Rdata')
raw_counts <- seurat__standard
save(raw_counts, file = 'site/counts_unfiltered.Rdata')

# extract counts and cpm
cpm <- RelativeCounts(integrated_obj@assays$RNA@counts, scale.factor= 1e6)


counts = integrated_obj@assays$RNA@counts
save(cpm, file = 'site/cpm.Rdata')
save(counts, file = 'site/counts.Rdata')

# extract meta filter
library(pool)
library(RSQLite)
db = args[6]
#db =  "site/scEiaD__2020_08_13__Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite"
pool <- dbPool(drv = SQLite(), dbname = db, idleTimeout = 3600000)
meta_filter <- pool %>% tbl('metadata_filter') %>% collect()
write_tsv(meta_filter, path = 'site/metadata_filter.tsv.gz')

# extract cell info
load('pipeline_data/cell_info/cell_info_labelled.Rdata')
write_tsv(cell_info_labels, path = 'site/cell_labels.tsv.gz')

# extract diff results
diff_results <- pool %>% tbl('PB_results') %>% collect()
readr::write_tsv(diff_results, file = 'pseudoBulk_diff_results.tsv.gz')
## write each PB_Test to separate file
## remove spaces in PB_Test as they will be used in file name export
diff_results$PB_Test <- gsub('\\s+','', diff_results$PB_Test)
diff_results %>% 
    group_by(PB_Test) %>%
    group_walk(~ readr::write_tsv(.x, paste0('site/pseudoBulkTable_', .y$PB_Test, ".tsv.gz")))
