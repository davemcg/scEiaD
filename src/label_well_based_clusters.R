
library(pool)
library(RSQLite)
library(tidyverse)
library(Seurat)


scEiaD <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/scEiaD__2020_08_20__Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)
load('/Volumes/data-1/projects/nei/mcgaughey/scEiaD_me/seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__sqrt__onlyWELL__batch__fastMNN__dims8__preFilter__mindist0.3__nneighbors30.umap.Rdata')
load('/Volumes/data-1/projects/nei/mcgaughey/scEiaD_me/umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__sqrt__onlyWELL__batch__fastMNN__dims8__preFilter__mindist0.3__nneighbors30.umap.Rdata')


meta_filter <- scEiaD %>% tbl('metadata_filter') %>% collect()


seurat_obj <- CreateSeuratObject(counts = integrated_obj@assays$RNA@counts)
for (i in colnames(umap)){
  seurat_obj <- AddMetaData(seurat_obj, umap[,i] %>% pull(1), col.name = i)
}

seurat_obj@reductions$mnn <- CreateDimReducObject(embeddings = integrated_obj@reductions$mnn@cell.embeddings,
                                                  loadings = integrated_obj@reductions$mnn@feature.loadings,
                                                  projected = integrated_obj@reductions$mnn@feature.loadings.projected,
                                                  key = 'mnn_',
                                                  assay = 'RNA')
seurat_obj@reductions$mnnUMAP <- CreateDimReducObject(embeddings = integrated_obj@reductions$mnnUMAP@cell.embeddings,
                                                      loadings = integrated_obj@reductions$mnnUMAP@feature.loadings,
                                                      projected = integrated_obj@reductions$mnnUMAP@feature.loadings.projected,
                                                      key = 'mnnUMAP_',
                                                      assay = 'RNA')


seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:8, reduction = 'mnn')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
#DimPlot(seurat_obj, reduction = 'mnnUMAP')

retina <- readr::read_tsv('https://raw.githubusercontent.com/davemcg/eyeMarkers/master/lists/retina_single_cell_markers__cowan2020.tsv')
library(tidyverse)
source(glue('{git_dir}src/make_gene_id_converter_table.R'))
#gene_id_converter <- read_tsv('~/data/massive_integrated_eye_scRNA/ensembl_biomart_human2mouse_macaque.tsv', skip = 1,
#                              col_names= c('hs_gene_id','hs_gene_id_v', 'mm_gene_id', 'mf_gene_id',
#                                           'hs_gene_name', 'mf_gene_name', 'mm_gene_name')) %>%
#  select(-hs_gene_id_v)

neb_plots <- list()
for (i in retina$CellType %>% unique()){
  print(i)
  genes <- retina %>% 
    filter(CellType == i) %>% 
    left_join(gene_id_converter %>% 
                mutate(Gene = hs_gene_name)) %>% 
    pull(hs_gene_id) %>% 
    unique() 
  genes <- genes[genes %in% row.names(seurat_obj)]
  try({neb_plots[[i]] <- plot_density(seurat_obj, genes, joint = TRUE, method = 'wkde')})
}


top_genes <- FindAllMarkers(seurat_obj,  group.by = 'seurat_clusters', logfc.threshold = 1, min.pct = 0.2)
top_genes_roc <- FindAllMarkers(seurat_obj,  group.by = 'seurat_clusters', logfc.threshold = 1, min.pct = 0.2, test.use = 'roc')

marker_exp_plots <- function(gene){
  if (grepl('ENSG', gene[1])){
    gene <- gene_id_converter %>% filter(hs_gene_id %in% (top_genes %>% filter(cluster == i, avg_logFC > 1) %>% pull(gene) %>% head(6))) %>% pull(hs_gene_name) %>% unique()
  }
  box_data <- scEiaD %>% tbl('grouped_stats') %>%
    filter(Gene %in% gene) %>%
    collect()
  grouping_features <- c('study_accession')
  exp_plot_facet <- 'CellType_predict'
  box_data <- box_data %>%
    #filter(!is.na(!!as.symbol(grouping_features))) %>%
    group_by_at(vars(one_of(c('Gene', exp_plot_facet, grouping_features)))) %>%
    summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
              cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
    full_join(., meta_filter %>%
                group_by_at(vars(one_of(grouping_features))) %>%
                summarise(Count = n())) %>%
    mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
    mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
           Expression = round(cpm * (`%` / 100), 2)) %>%
    select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
    arrange(-Expression) %>%
    rename(`Cells # Detected` = cell_exp_ct,
           `Total Cells` = Count,
           `Mean CPM` = Expression,
           `% of Cells Detected` = `%`) %>%
    tidyr::drop_na()
  box_data$Group <- box_data[,c(2:(length(grouping_features)+1))] %>% tidyr::unite(x, sep = ' ') %>% pull(1)
  
  box_data %>%
    filter(!grepl('Hufnagel', study_accession)) %>% 
    ggplot(aes(x=Gene, y = `Mean CPM`, color = !!as.symbol(grouping_features))) +
    geom_boxplot(color = 'black', outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(aes(size = `Total Cells`), grouponX = TRUE) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()), 20)) +
    theme(legend.position="bottom") +
    facet_wrap(ncol = 15, scales = 'free_x', vars(!!as.symbol(exp_plot_facet))) +
    ggtitle(paste0('Cluster ', i))
  
}

marker_diff_test <- function(gene){
  if (grepl('ENSG', gene[1])){
    gene <- gene_id_converter %>% filter(hs_gene_id %in% (top_genes %>% filter(cluster == i, avg_logFC > 1) %>% pull(gene) %>% head(6))) %>% pull(hs_gene_name) %>% unique()
  }
  results <- scEiaD %>% tbl("PB_results") %>% 
    filter(Gene %in% gene) %>% collect() %>% 
    filter(grepl('CellType \\(Pre', PB_Test)) %>% 
    arrange(FDR) %>% 
    head(50)
  results
}


exp_plots <- list()
marker_diff <- list()
marker_auc <- list()
for (i in top_genes$cluster %>% unique()){
  print(i)
  exp_plots[[i]] <- marker_exp_plots(top_genes %>% filter(cluster == i, avg_logFC > 1) %>% pull(gene) %>% head(6))
  marker_diff[[i]] <- marker_diff_test(top_genes %>% filter(cluster == i, avg_logFC > 1) %>% pull(gene) %>% head(6))
  marker_auc[[i]] <- marker_diff_test(top_genes_roc %>% filter(cluster == i, avg_logFC > 1) %>% pull(gene) %>% head(6))
  ggsave(plot = exp_plots[[i]], filename = paste0('plots/exp_plot_cluster_', i, '.png'), device = 'png', width = 20, height = 20, dpi = 'retina')
}


ct2 <- seurat_obj@meta.data %>% as_tibble() %>% mutate(CellType = case_when(seurat_clusters == '0' ~ 'Rods',
                                                                            seurat_clusters == '1' ~ 'Muller Glia',
                                                                            seurat_clusters == '2' ~ 'RPCs',
                                                                            seurat_clusters == '3' ~ 'RPCs',
                                                                            seurat_clusters == '4' ~ 'Cones',
                                                                            seurat_clusters == '5' ~ 'Bipolar Cells',
                                                                            seurat_clusters == '6' ~ 'RPC',
                                                                            seurat_clusters == '7' ~ 'Amacrine Cells',
                                                                            seurat_clusters == '8' ~ 'Rods',
                                                                            seurat_clusters == '9' ~ 'Mesenchymal/RPE/Endothelial',
                                                                            seurat_clusters == '10' ~ 'Bipolar Cells',
                                                                            seurat_clusters == '11' ~ 'Amacrine Cells',
                                                                            seurat_clusters == '12' ~ 'Bipolar Cells',
                                                                            seurat_clusters == '13' ~ 'Rods',
                                                                            TRUE ~ 'RPCs')) %>% 
  pull(CellType)



seurat_obj <- AddMetaData(seurat_obj, ct2, col.name = 'CellType')
seurat_obj <- AddMetaData(seurat_obj, ct2, col.name = 'CellType2')
DimPlot(seurat_obj, reduction = 'mnnUMAP', group.by = 'CellType')



save(seurat_obj, file = 'data/well_data_seurat_obj_labelled.Rdata')

