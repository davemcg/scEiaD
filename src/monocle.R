# expression_matrix, a numeric matrix of expression values, where rows are genes, and columns are cells
# cell_metadata, a data frame, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
# gene_metadata, an data frame, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

# The expression value matrix must:
#   have the same number of columns as the cell_metadata has rows.
#   have the same number of rows as the gene_metadata has rows.

# Additionally:
#   row names of the cell_metadata object should match the column names of the expression matrix.
#   row names of the gene_metadata object should match row names of the expression matrix.
#   one of the columns of the gene_metadata should be named "gene_short_name", 
#     which represents the gene symbol or simple name (generally used for plotting) for each gene.
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(tidyverse)
library(monocle3)
#load('seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims200__mindist0.1__nneighbors50.umap.Rdata')
load(args[1])
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__counts__full__batch__scVI__dims200__mindist0.1__nneighbors50.umap.Rdata')
load(args[2])


plot_file <- args[3]
plot_file_no_label <- args[4]
output <- args[5]

# gene_annotation <- grep('^MT-', integrated_obj@assays$RNA@var.features, invert =TRUE, value = TRUE) %>% 
#   enframe(name = NULL, value ='gene_short_name') %>% 
#   data.frame()
# row.names(gene_annotation) <- gene_annotation$gene_short_name
gene_annotation <- row.names(integrated_obj) %>% 
  enframe(name = NULL, value ='gene_short_name') %>% 
  data.frame()
row.names(gene_annotation) <- row.names(integrated_obj)

cell_metadata <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(umap %>% select(Barcode, study_accession:Method), by = 'Barcode') %>% 
  data.frame()
row.names(cell_metadata) <- cell_metadata$Barcode

counts <- integrated_obj@assays$RNA@counts
#row.names(counts) <- gene_annotation$gene_short_name
#colnames(counts) <- row.names(integrated_obj@meta.data)



cds_retina <- new_cell_data_set(counts,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# cds_retina@clusters$UMAP$partitions <- meta$cluster_knn10
# names(cds_retina@clusters$UMAP$partitions) <- row.names(integrated_obj@meta.data)

cds_retina@int_colData$reducedDims$PCA <- as.matrix(integrated_obj@reductions[["scVI"]]@cell.embeddings)
cds_retina@int_colData$reducedDims$UMAP <- as.matrix(integrated_obj@reductions[["scviUMAP"]]@cell.embeddings)

# monocle uses a igrpah structure with some more data I haven't computed?
# for now re-run
cds_retina <- cluster_cells(cds_retina, method = 'PCA')
cds_retina <- learn_graph(cds_retina)

rownames(cds_retina@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds_retina@int_colData$reducedDims$UMAP) <- NULL

png(filename = plot_file, width = 5000, height = 4000, res = 400)
plot_cells(cds_retina, color_cells_by = 'CellType_predict', 
           label_roots=FALSE, 
           label_leaves = FALSE, 
           alpha = 0.8, 
           label_branch_points = FALSE, 
           group_label_size= 5)
dev.off()

png(filename = plot_file_no_label, width = 9000, height = 4000, res = 400)
plot_cells(cds_retina, color_cells_by = 'CellType_predict', 
           label_roots=FALSE, 
           label_leaves = FALSE, 
           alpha = 0.8, 
           label_branch_points = FALSE, 
           group_label_size= 5,
           label_cell_groups = FALSE)
dev.off()

save(cds_retina, 
     file = output)
# 
# cds_retina@reduce_dim_aux[["UMAP"]] = as.matrix(integrated_obj@reductions[["scviUMAP"]]@cell.embeddings)
# 
# cds_retina@preprocess_aux$gene_loadings <- integrated_obj@reductions$scVI@cell.embeddings

