library(Matrix)
library(loomR)

args <- commandArgs(trailingOnly = TRUE)
load(args[1])

create(filename= args[2],
           overwrite = TRUE,
           data = matrix,
           cell.attrs = list(batch = loom_batch,
                             batch_indices = loom_batch_indices))

