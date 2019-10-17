# cluster purity plot

# the idea is that, theoretically..., every seurat cluster is one cell type
# so let's quantify that with the known labels
args <- commandArgs(trailingOnly = TRUE)

library(cowplot)
library(tidyverse)
library(splancs)
library(ggsci)

umap_folder <- args[1]
umap_files <- list.files(umap_folder, "*downsample*", full.names = TRUE)
all_obj <- list()
for (i in umap_files){
  load(i)
  set <- gsub('umap/|.umap.Rdata','', i)
  all_obj[[set]] <- umap
}

output_file1 <- args[2]
output_file2 <- args[3]
output_file3 <- args[4]

cluster_purity_plot <- function(obj, algorithm_name, 
                                grouping_var_1 = 'seurat_clusters', 
                                grouping_var_2 = 'CellType'){
  obj <- obj %>% filter(!is.na(!!(sym(grouping_var_2))))
  p1_mean <- obj %>% group_by(!!(sym(grouping_var_1)), !!(sym(grouping_var_2))) %>% 
    summarise(Count = n()) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(!!(sym(grouping_var_1))) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
    pull(Count) %>% mean() %>% round(., 2)
  p1_cluster_count <- obj %>% group_by(!!(sym(grouping_var_1)), !!(sym(grouping_var_2))) %>% 
    summarise(Count = n()) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(!!(sym(grouping_var_1))) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) %>% 
    pull(!!(sym(grouping_var_1))) %>% unique() %>% length()
  obj <- obj %>% group_by(!!(sym(grouping_var_1)), !!(sym(grouping_var_2))) %>% 
    summarise(Count = n(),  Grouping2 = list(unique(!!(sym(grouping_var_2))))) %>% 
    mutate(Perc = Count/sum(Count)) %>% 
    filter(Perc > 0.1) %>% 
    ungroup() %>% group_by(!!(sym(grouping_var_1))) %>% 
    summarise(Total = sum(Count), Count = n()) %>% filter(Total > 500) 
  obj$Info <- algorithm_name
  obj <- obj %>% 
    separate(Info, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', 'Dims'), sep = '__')
  cells_per_group <- obj %>% 
    group_by(Count) %>% 
    summarise(Sum = sum(Total)) %>% 
    mutate(Count = case_when(Count >=2 ~ '2 or more',
                             TRUE ~ '1')) %>% 
    ungroup() %>% 
    mutate(Total = sum(Sum)) %>% 
    group_by(Count) %>% 
    summarise(Ratio = sum(Sum) / max(Total))
  list(mean = p1_mean, cluster_count = p1_cluster_count, dist = obj, cells_per_group = cells_per_group)
}

# only keep full 
full <- grep('full', names(all_obj), value = TRUE)
full_obj <- all_obj[full]
bipolar_muller_rod <- full_obj %>% map(function(x) mutate(x, CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', CellType) )) %>% 
  map(function(x) filter(x, CellType %in% c('Bipolar Cells', 'Muller Glia', 'Rods'))) %>% 
  map(function(x) filter(x, Age  > 10)) %>% 
  map(function(x) mutate(x, celltype_study = paste0(CellType, '__', study_accession)))


purity_calcs_cluster_vs_celltype <- bipolar_muller_rod %>% imap(cluster_purity_plot, "seurat_clusters", "CellType")
purity_calcs_cluster_vs_study_and_cell_type <- bipolar_muller_rod %>% imap(cluster_purity_plot, "seurat_clusters", "celltype_study")
#purity_calcs_cluster_vs_age <- bipolar_muller_rod %>% imap(cluster_purity_plot, "seurat_clusters", "Age")


dist_plotter <- function(output_file, obj, title = "Distribution of Purity of Clusters by Cell Type"){
  pdf(output_file, width = 20, height = 4)
  print(obj %>% 
          map(3) %>% # third item in list is the purity distributions
          bind_rows() %>% 
          ggplot(aes(x=Method, y = Count, colour = Transform)) + 
          geom_violin(draw_quantiles = c(0.5)) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          facet_wrap(~Covariate + Set + Dims, nrow = 1) +
          ggtitle(title) 
  )
  dev.off()   
  #  }
  
}

mean_plotter <- function(output_file, obj, title = "Mean Purity of Clusters by Cell Type"){
  pdf(output_file, width = 10, height = 8)
  print(cbind(obj %>% 
                map(1) %>% # first item is the mean
                enframe(value = 'Mean') %>% 
                unnest(Mean),
              obj %>% 
                map(2) %>% # second is the number of clusters
                enframe(value = 'Count') %>% unnest(Count) %>% select(Count)) %>% 
          separate(name, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', "Dims"), sep = '__') %>% 
          as_tibble() %>% 
          ggplot(aes(x=Method, y = Mean, colour = Transform)) + 
          geom_point() + 
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          facet_wrap(~CellType, nrow = 2) +
          ggsci::scale_color_aaas() +
          ggtitle(title))
  dev.off()
  #}
}

ratio_plotter <- function(output_file, obj, title = "Purity of Clusters by Cell Type"){
  pdf(output_file, width = 10, height = 8)
  print(obj %>% 
          map(4) %>% # second is the number of clusters
          bind_rows(.id = 'name') %>% 
          separate(name, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', "Dims"), sep = '__') %>% 
          as_tibble() %>% 
          ggplot(aes(x=Method, y = Ratio, fill = Count)) + 
          geom_bar(stat = 'identity') +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          facet_grid(vars(Method), vars(Dims)) +
          ggsci::scale_color_aaas() +
          ggtitle(title))
  dev.off()
  #}
}

extract <- function(obj1, obj2){
  rbind(obj1 %>% 
          map(4) %>% # second is the number of clusters
          bind_rows(.id = 'name') %>% 
          separate(name, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', "Dims"), sep = '__') %>% 
          as_tibble() %>% arrange(Ratio) %>% filter(Count == 1) %>% rename(Score = Ratio) %>% mutate(Type = 'CellType Purity'),
        obj2 %>% map(4) %>% # second is the number of clusters
          bind_rows(.id = 'name') %>% 
          separate(name, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', "Dims"), sep = '__') %>% 
          as_tibble() %>% arrange(Ratio) %>% filter(Count != 1) %>% rename(Score = Ratio) %>% mutate(Type = 'Study Blending'))
}

extract(purity_calcs_cluster_vs_celltype, 
        purity_calcs_cluster_vs_study_and_cell_type) %>% 
  as_tibble() %>% 
  mutate(Method = factor(Method, levels = c('fastMNN', 'CCA', 'harmony', 'scanorama', 'combat', 'none')),
         Type = factor(Type, levels = c('Study Blending', 'CellType Purity')),
         Dims = gsub('dims','', Dims) %>% as.numeric()) %>% 
  ggplot(aes(x=Method, y = Score, fill = as.factor(Dims), color = Type)) + 
  geom_bar(stat = 'identity', size = 2) + 
  facet_grid(vars(Covariate), vars(Dims)) + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggsci::scale_fill_aaas() + scale_colour_manual(values = c("gray","yellow"))

# draw a convex hull around each seurat_cluster and plot size
# ideally want each cluster to be "tight" (i.e. not a big ol blob)
cluster_area <- function(obj, cluster_col = 'seurat_clusters'){
  obj <- obj[,c('UMAP_1', 'UMAP_2', cluster_col)]
  obj <- obj %>% mutate(col = as.numeric(as.character(!!sym(cluster_col))))
  
  areas <- c()
  for (i in obj %>% pull(col) %>% unique()){
    f_obj <- obj %>% filter(col == i)
    # skip clustesr with < 500 cells
    if (nrow(f_obj) >= 500) {
      # convex hull
      outer_coords <- chull(f_obj$UMAP_1, f_obj$UMAP_2)
      # Calculate area in convex hull
      area <- splancs::areapl(f_obj[outer_coords,c('UMAP_1', 'UMAP_2')] %>% as.matrix()) 
      areas <- c(areas, area)
    }
  }
  areas %>% enframe()
}


#dist_plotter('test1', purity_calcs_cluster_vs_celltype, "Distribution of Purity of Clusters by Cell Type")
#dist_plotter('test2', purity_calcs_cluster_vs_study_and_cell_type, "Distribution of Purity of Clusters by Study")
#dist_plotter('test3', purity_calcs_cluster_vs_age, "Distribution of Purity of Clusters by Age")

mean_plotter(args[2], purity_calcs_cluster_vs_celltype, "Mean Purity of Clusters by Cell Type (Ideal is 1)")
mean_plotter(args[3], purity_calcs_cluster_vs_study_and_cell_type, "Diversity of Cell Type Clusters by Study (Higher is Better)")
#mean_plotter('test6', purity_calcs_cluster_vs_age, "Mean Purity of Clusters by Age")

pdf(args[4], width = 6, height = 6)
area_plot <- full_obj %>% map(cluster_area) %>% bind_rows(.id='stuff') %>% separate(stuff, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', "Dims"), sep = '__')
area_plot %>% 
  group_by(Transform, Covariate, Method, Dims) %>% 
  summarise(value = median(value)) %>% 
  ggplot(aes(x=Method, y = value, colour = Transform)) + 
  geom_point(stat = 'identity') + 
  theme_minimal() + 
  ggsci::scale_color_aaas() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~Covariate + Dims, nrow = 2) + 
  ylab('Median Cluster Area') + ggtitle("fastMNN has the smallest median cluster size")
dev.off()