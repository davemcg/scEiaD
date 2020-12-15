library(tidyverse)
args <- commandArgs(trailingOnly = T)
save(args, file='testing/mg_args.Rdata')
mouse_gtf <- rtracklayer::readGFF(args[1])
maca_gtf <- rtracklayer::readGFF(args[2])
human_gtf <- rtracklayer::readGFF(args[3])

mouse_outfile <- args[5]
maca_outfile <- args[6]
human_outfile <- args[7]
human_mito_genes <- human_gtf %>% filter(type == 'gene', grepl('^MT-', gene_name)) %>% pull(gene_id) %>% paste0('.')
mouse_mito_genes <- mouse_gtf %>% filter(type == 'gene', grepl('^mt-', gene_name)) %>% pull(gene_id) %>% paste0('.')
maca_mito_genes <- maca_gtf %>% filter(type == 'transcript', seqid == 'MT', !is.na(transcript_name)) %>% pull(gene_id) %>% unique
#weirdly, the human and mouse names get '.' added, but maca does not fo

write(human_mito_genes, file = human_outfile, sep = '\n')
write(mouse_mito_genes, file = mouse_outfile,sep= '\n')
write(maca_mito_genes, file = maca_outfile,sep= '\n')
