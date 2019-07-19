#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(purrr)
library(readr)
library(stringr)
library(dplyr)
library(data.table)

files <- args[seq(3,length(args))]
SRS = str_extract(files, 'SRS.*/') %>% gsub('/','',.)
tpm <- files %>% 
  map(fread, select = c('target_id', 'tpm')) %>% 
  reduce(left_join, by = 'target_id')
count <- files %>% 
  map(fread, select = c('target_id', 'est_counts')) %>% 
  reduce(left_join, by = 'target_id')

target_lengths <- fread(files[1], select = c('target_id','length','eff_length'))

tpm <- left_join(target_lengths, tpm, by = 'target_id')
count <- left_join(target_lengths, tpm, by = 'target_id') 

colnames(tpm) <- c('target_id', 'length', 'eff_length', SRS)
colnames(count) <- c('target_id', 'length', 'eff_length', SRS)

sparse_tpm <- tpm[,2:ncol(tpm)]%>% as.matrix() %>% Matrix::Matrix(., sparse = TRUE)
row.names(sparse_tpm) <- tpm[,1]
sparse_count <- count[,2:ncol(count)]%>% as.matrix() %>% Matrix::Matrix(., sparse = TRUE)
row.names(sparse_count) <- count[,1]


save(sparse_tpm, file = args[1])
save(sparse_count, file = args[2])
