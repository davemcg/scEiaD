#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tximport)
library(readr)
library(stringr)
library(dplyr)
library(Matrix)

files <- args[seq(4,length(args))]

tx2gene <- read_tsv(args[3], col_names = F) %>% select(X1, X3)

SRS = str_extract(files, 'SRS.*/') %>% gsub('/','',.)

txi_count <- tximport(files, type = 'kallisto', tx2gene = tx2gene %>% select(X1, X3), countsFromAbundance = 'no')

count <- txi_count$counts
colnames(count) <- SRS
count <- count %>% Matrix(., sparse = TRUE)

save(count, file = args[1])

txi_count <- tximport(files, type = 'kallisto', txOut = TRUE, countsFromAbundance = 'no')

count <- txi_count$counts
colnames(count) <- SRS
count <- count %>% Matrix(., sparse = TRUE)
save(count, file = args[2])
