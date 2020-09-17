library(tidyverse)
library(ggraph)
library(tidygraph)

args = commandArgs(trailingOnly=TRUE)

load(args[1])
out <- args[2]

centrality1 <- tidygraph::as_tbl_graph(mst_celltype[[1]]) %>% 
    activate('nodes') %>% 
  mutate(type = case_when(grepl('RPC|Neuro|Precu', name) ~ 'Progenitor',
                          grepl('Rods|Cones', name) ~ 'PR',
                          grepl('Chorio|Endothe|Smooth|Peric', name) ~ 'Vaculature',
                          grepl('glia|Astro|Schw', name) ~ 'Glia',
                          grepl('ipolar', name) ~ 'Bipolar Cells',
                          grepl('Muller', name) ~ 'Muller Glia',
                          grepl('Amacrine', name) ~ 'Amacrine',
                          grepl('Retinal Ganglion', name) ~ 'Retinal Ganglion',
                          grepl('Melano', name) ~ 'Melanocytes',
                          grepl('Horizon', name) ~ 'Horizontal',
                          TRUE ~ name),
         ac = centrality_authority(),
         bc = centrality_betweenness(),
         cc = centrality_closeness(),
         assortativity = graph_assortativity(type, directed = FALSE)) %>% 
	as_tibble() %>% 
	mutate(Group = 'CellType')

centrality2 <- tidygraph::as_tbl_graph(mst_cluster[[1]]) %>%
    activate('nodes') %>%
      mutate(type = case_when(grepl('RPC|Neuro|Precu', name) ~ 'Progenitor',
                          grepl('Rods|Cones', name) ~ 'PR',
                          grepl('Chorio|Endothe|Smooth|Peric', name) ~ 'Vaculature',
                          grepl('glia|Astro|Schw', name) ~ 'Glia',
                          grepl('ipolar', name) ~ 'Bipolar Cells',
                          grepl('Muller', name) ~ 'Muller Glia',
                          grepl('Amacrine', name) ~ 'Amacrine',
                          grepl('Retinal Ganglion', name) ~ 'Retinal Ganglion',
                          grepl('Melano', name) ~ 'Melanocytes',
                          grepl('Horizon', name) ~ 'Horizontal',
                          TRUE ~ name),
         ac = centrality_authority(),
         bc = centrality_betweenness(),
         cc = centrality_closeness(),
         assortativity = graph_assortativity(type, directed = FALSE)) %>%
	as_tibble() %>% 
	mutate(Group = 'Cluster')

centrality <- bind_rows(centrality1, centrality2)
write_tsv(centrality, path = out)

