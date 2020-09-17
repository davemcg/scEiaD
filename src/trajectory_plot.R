library(tidyverse)
library(cowplot)
library(ggraph)
library(tidygraph)

args = commandArgs(trailingOnly=TRUE)
load(args[1])
out <- args[2]
graph_celltype <- tidygraph::as_tbl_graph(mst_celltype[[1]]) %>% 
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
			bc = centrality_betweenness()) %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link() + 
  geom_node_point(aes(colour = type, size = bc)) +
  geom_node_text(aes(label = name), colour = 'black', vjust = 0.4) +
  theme_graph() +
  scale_color_manual(values = pals::polychrome() %>% unname())

graph_cluster <- tidygraph::as_tbl_graph(mst_cluster[[1]]) %>%
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
		bc = centrality_betweenness()) %>% 
  ggraph(layout = 'kk') +
  geom_edge_link() +
  geom_node_point(aes(colour = type, size = bc)) +
  geom_node_text(aes(label = name), colour = 'black', vjust = 0.4) +
  theme_graph() +
  scale_color_manual(values = pals::polychrome() %>% unname())


png(out, width = 2000, height = 4000, res = 150)
plot_grid(graph_celltype, graph_cluster, ncol = 1)
dev.off()


