#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(purrr)
library(readr)
library(stringr)
library(dplyr)

files <- c('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/quant/SRS3106638/abundance.tsv.gz',
           '/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/quant/SRS3106639/abundance.tsv.gz',
           '/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/quant/SRS3106675/abundance.tsv.gz',
           '/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/quant/SRS3106675/abundance.tsv.gz')

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

save(tpm, file = args[1])
save(count, file = args[2])