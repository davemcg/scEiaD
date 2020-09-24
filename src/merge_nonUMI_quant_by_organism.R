#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
save(args, file = 'testing/spnew_mnumi_arg.Rdata')
library(tximport)
library(readr)
library(stringr)
library(dplyr)
library(Matrix)
files <- args[seq(5,length(args))]
system(paste0('mkdir -p ', dirname(args[1])))

tx2gene <- read_tsv(args[3], col_names = c('transcript_id_dot', 'gene_id_dot')) %>% 
  mutate(gene_id = str_remove(gene_id_dot, '\\.$'), 
         transcript_id = str_remove(transcript_id_dot, '\\.$') )
gtf <- rtracklayer::readGFF(args[4]) %>% as_tibble 
geneid2gene_name <- gtf %>% filter(type == "transcript") %>% select(gene_id, gene_name) %>% distinct
tx2gene <- tx2gene %>% 
  left_join( geneid2gene_name) %>% 
  mutate(gene_name = replace(gene_name, is.na(gene_name), gene_id[is.na(gene_name)]))
#save(tx2gene, file = args[3])

SRS = str_split(files, 'quant/') %>% sapply(function(x) x[2]) %>% str_split('/') %>% sapply(function(x) x[1])

txi_count <- tximport(files, type = 'kallisto', tx2gene = tx2gene %>% select(transcript_id_dot, gene_id), countsFromAbundance = 'no')
count <- txi_count$counts
colnames(count) <- SRS
rownames(count) <- str_remove(rownames(count), '\\.$') 
count <- count %>% Matrix(., sparse = TRUE)

save(count, file = args[1])

txi_count <- tximport(files, type = 'kallisto', txOut = TRUE, countsFromAbundance = 'no')

count <- txi_count$counts
colnames(count) <- SRS
count <- count %>% Matrix(., sparse = TRUE)
save(count, file = args[2])
