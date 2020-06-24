library(monocle3)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
threads = as.numeric(args[2])
comp = args[3]
out = args[4]

topN <- top_markers(cds_retina, group_cells_by= comp, reference_cells=150000, cores = threads, verbose = FALSE)

save(topN, file = out)
