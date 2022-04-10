load('seurat_obj/integrated/n_features-6000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-10__epochs-15__dist-0.1__neighbors-50.umap.Rdata')

int_obj_A <- integrated_obj

load('seurat_obj/integrated/n_features-6000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-10__epochs-15__preFilter.seuratV3.Rdata')


load('pipeline_data/clean_quant/all_species_full_sparse_matrix_intron.Rdata')

cells <- colnames(integrated_obj)
genes <- row.names(integrated_obj)
intron_mat <- all_cells_all_species_matrix_intron[genes, cells]
#integrated_obj[['spliced']] <- integrated_obj[['RNA']]
integrated_obj[['unspliced']] <- CreateAssayObject(counts = intron_mat)
# find var features for unspliced
integrated_obj@active.assay <- 'unspliced'
integrated_obj <- FindVariableFeatures(integrated_obj, selection.method = 'vst')
# reset back to RNA as active assay (h5ad will use this as "X", so can't leave unspliced as default)
integrated_obj@active.assay <- 'RNA'
# integrated_obj[['unspliced']]@var.features <- integrated_obj[['RNA']]@var.features
# # add feature selection genes info to unspliced
# var_f <- var_f[,'vst.variable', drop = FALSE]
# integrated_obj@assays$unspliced@var.features <- var_f
# copy over umap 2D
integrated_obj[['umap']] <- int_obj_A[['scviUMAP']]
integrated_obj@reductions$umap@key <- 'umap_'

system('mkdir -p pipeline_data/velocity')
save(integrated_obj, file = 'pipeline_data/velocity/n_features-6000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-10__epochs-15__preFilter.addIntron.seuratV3.Rdata', compress = FALSE)
