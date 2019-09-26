# cluster purity plot

# the idea is that, theoretically..., every seurat cluster is one cell type
# so let's quantify that with the known labels
args <- commandArgs(trailingOnly = TRUE)

library(cowplot)
library(tidyverse)

umap_folder <- args[1]
umap_files <- list.files(umap_folder, "*umap.Rdata", full.names = TRUE)
all_obj <- list()
for (i in umap_files){
  load(i)
  set <- gsub('umap/|.umap.Rdata','', i)
  all_obj[[set]] <- umap
}

output_file1 <- args[2]
output_file2 <- args[3]

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
  list(mean = p1_mean, cluster_count = p1_cluster_count, dist = obj)
}

purity_calcs_cluster_vs_celltype <- all_obj %>% imap(cluster_purity_plot, "seurat_clusters", "CellType")
purity_calcs_cluster_vs_study <- all_obj %>% imap(cluster_purity_plot, "seurat_clusters", "study_accession")
purity_calcs_cluster_vs_age <- all_obj %>% imap(cluster_purity_plot, "seurat_clusters", "Age")


dist_plotter <- function(output_prefix, obj, title = "Distribution of Purity of Clusters by Cell Type"){
  for (i in c("early", "late", "full")){
    subset = grep(i, names(obj), value = TRUE)
    pdf(paste0(output_prefix, '_', i, '.pdf'), width = 20, height = 4)
    #pdf('test.pdf', width = 20, height = 4)
    print(obj[subset] %>% 
            map(3) %>% # third item in list is the purity distributions
            bind_rows() %>% 
            ggplot(aes(x=Method, y = Count, colour = Transform)) + 
            geom_violin(draw_quantiles = c(0.5)) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            facet_wrap(~Covariate + Set + Dims, nrow = 1) +
            ggtitle(title) 
    )
    dev.off()   
  }
  
}

mean_plotter <- function(output_prefix, obj, title = "Mean Purity of Clusters by Cell Type"){
  for (i in c("early", "late", "full")){
    subset = grep(i, names(obj), value = TRUE)
    pdf(paste0(output_prefix, '_', i, '.pdf'), width = 20, height = 4)
    print(cbind(obj[subset] %>% 
                  map(1) %>% # first item is the mean
                  enframe(value = 'Mean') %>% 
                  unnest(Mean),
                obj[subset] %>% 
                  map(2) %>% # second is the number of clusters
                  enframe(value = 'Count') %>% unnest(Count) %>% select(Count)) %>% 
            separate(name, into = c('Species', 'Transform', 'Set', 'Covariate', 'Method', "Dims"), sep = '__') %>% 
            as_tibble() %>% 
            ggplot(aes(x=Method, y = Mean, colour = Transform, size = Count)) + 
            geom_point() + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            facet_wrap(~Covariate + Set + Dims, nrow = 1) +
            ggtitle(title))
    dev.off()
  }
}


dist_plotter('test1', purity_calcs_cluster_vs_celltype, "Distribution of Purity of Clusters by Cell Type")
dist_plotter('test2', purity_calcs_cluster_vs_study, "Distribution of Purity of Clusters by Study")
dist_plotter('test3', purity_calcs_cluster_vs_age, "Distribution of Purity of Clusters by Age")

mean_plotter('test4', purity_calcs_cluster_vs_celltype, "Mean Purity of Clusters by Cell Type")
mean_plotter('test5', purity_calcs_cluster_vs_study, "Mean Purity of Clusters by Study")
mean_plotter('test6', purity_calcs_cluster_vs_age, "Mean Purity of Clusters by Age")

