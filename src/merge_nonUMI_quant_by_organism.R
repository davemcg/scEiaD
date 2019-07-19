#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tximport)
library(readr)
library(stringr)
library(dplyr)


files <- args[seq(3,length(args))]

tx2gene <- read_tsv(args[2], col_names = F) %>% select(X1, X3)

SRS = str_extract(files, 'SRS.*/') %>% gsub('/','',.)

txi_tpm <- tximport(files, type = 'kallisto', tx2gene = tx2gene %>% select(X1, X3), countsFromAbundance = 'lengthScaledTPM')

tpm <- txi_tpm$counts
colnames(tpm) <- SRS
tpm <- tpm %>% Matrix(., sparse = TRUE)

save(tpm, file = args[1])