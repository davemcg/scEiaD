args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(Matrix)

load(args[1]) # matrix
mat_name <- args[2] # get name of matrix

mat <- get(mat_name)

drop_samples <- colnames(mat) %>% str_extract(., '_\\w+') %>% gsub('_', '', .)
drop_samples <- drop_samples[!is.na(drop_samples)]

well_samples <- grep('_', colnames(mat), value=TRUE, invert = TRUE)

system('mkdir -p site/sample_matrices')
matrix_writer <- function(big_matrix, sample){
    cells <- grep(sample, colnames(big_matrix), value = TRUE)
    sample_matrix <- big_matrix[,cells]
    save(sample_matrix, file = paste0('site/sample_matrices/', sample, '.Rdata'))
}

for (i in unique(drop_samples)){
    print(i)
    matrix_writer(mat, i)
}

org <- str_split(mat_name, '_')[[1]][1]
if (length(well_samples) > 0){
	sample_matrix <- mat[,well_samples]
	save(sample_matrix, file = paste0('site/sample_matrices/', org, '_well.Rdata'))
}
