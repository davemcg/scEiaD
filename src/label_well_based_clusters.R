library(pool)
library(RSQLite)
library(tidyverse)
library(Seurat)


scEiaD <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)

meta_filter <- scEiaD %>% tbl('metadata_filter')

load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__libSize__onlyWELL__batch__fastMNN__dims30__preFilter__mindist0.1__nneighbors50.seurat.Rdata')
load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features2000__libSize__onlyWELL__batch__fastMNN__dims30__preFilter__mindist0.1__nneighbors50.umap.Rdata')




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
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = 'mnn')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
#DimPlot(seurat_obj, reduction = 'mnnUMAP')


DimPlot(seurat_obj, reduction = 'mnnUMAP', group.by = 'seurat_clusters')

# seurat_obj <- SetIdent(object = seurat_obj,value =  "cluster")

top_genes <- FindAllMarkers(seurat_obj,  group.by = 'seurat_clusters', logfc.threshold = 1, min.pct = 0.2)
top_genes_roc <- FindAllMarkers(seurat_obj,  group.by = 'seurat_clusters', logfc.threshold = 1, min.pct = 0.2, test.use = 'roc')


meta_filter <- scEiaD_2020_v01 %>% tbl("metadata_filter") %>% collect()
marker_exp_plots <- function(gene){
  box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
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
  results <- scEiaD_2020_v01 %>% tbl("PB_results") %>% 
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
  ggsave(plot = exp_plots[[i]], filename = paste0('plots/exp_plot_cluster_', i, '.png'), device = 'png', width = 20, dpi = 'retina')
}


ct <- seurat_obj@meta.data %>% as_tibble() %>% mutate(CellType = case_when(seurat_clusters == '0' ~ 'Cones',
                                                                           seurat_clusters == '1' ~ 'Rods',
                                                                           seurat_clusters == '2' ~ 'Bipolar Cells',
                                                                           seurat_clusters == '3' ~ 'Retinal Ganglion Cells',
                                                                           seurat_clusters == '4' ~ 'Rods',
                                                                           seurat_clusters == '5' ~ 'RPCs',
                                                                           seurat_clusters == '6' ~ 'RPCs',
                                                                           seurat_clusters == '7' ~ 'Muller Glia',
                                                                           seurat_clusters == '8' ~ 'Muller Glia',
                                                                           seurat_clusters == '9' ~ 'Muller Glia Progenitor',
                                                                           seurat_clusters == '10' ~ 'Muller Glia',
                                                                           seurat_clusters == '11' ~ 'Bipolar Cells',
                                                                           seurat_clusters == '12' ~ 'Bipolar Cells',
                                                                           seurat_clusters == '13' ~ 'Amacrine Cells',
                                                                           seurat_clusters == '14' ~ 'Amacrine Cells')) %>% 
  pull(CellType)

seurat_obj <- AddMetaData(seurat_obj, ct, col.name = 'CellType')
DimPlot(seurat_obj, reduction = 'mnnUMAP', group.by = 'CellType')

save(seurat_obj, file = 'data/well_data_seurat_obj_labelled.Rdata')
