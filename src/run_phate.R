library(data.table)
library(phateR)
library(tidyverse)

print(reticulate::py_discover_config("phate"))

args <- commandArgs(trailingOnly = TRUE)

load(args[1])

phate_2D <- phate(normalized_values, knn = 100, n.jobs = 24)
save(phate_2D, file = args[2])
