library(tidyverse)

files_drop <- list.files('quant/', pattern = '^matrix.seu*', recursive=TRUE, full.names=TRUE)
files_well <- list.files('quant/', pattern = 'abundance_spliced.tsv.gz', recursive=TRUE, full.names=TRUE)

drop_meta <- list()
for (i in files_drop){
	load(i)
	drop_meta[[i]] <- seu@meta.data %>% as_tibble(rownames = 'Barcode')
	drop_meta[[i]]$sample <- str_extract(i, 'SRS\\d+|ERS\\d+|E-MTAB\\d+|iPSC_RPE_scRNA_\\d+')
}

# grab mito genes
gtf_HS <-rtracklayer::readGFF('references/gtf/hs-homo_sapiens_anno.gtf.gz')
gtf_MM <-rtracklayer::readGFF('references/gtf/mm-mus_musculus_anno.gtf.gz')
gtf_MF <-rtracklayer::readGFF('references/gtf/mf-macaca_mulatta_anno.gtf.gz')

mito_genes_HS <- gtf_HS %>% filter(type == 'transcript', grepl('^MT-', transcript_name, ignore.case = TRUE)) %>% pull(transcript_id) %>% paste0('.')
mito_genes_MM <- gtf_MM %>% filter(type == 'transcript', grepl('^MT-', transcript_name, ignore.case = TRUE)) %>% pull(transcript_id) %>% paste0('.')
mito_genes_MF <-   gtf_MF %>% filter(type == 'transcript', seqid == 'MT', !is.na(transcript_name)) %>% pull(transcript_id)

mito_genes <- c(mito_genes_HS, mito_genes_MM, mito_genes_MF)



well_meta <- list()
inc <- 1
for (i in files_well){
	print(inc); inc = inc + 1
	sample <- str_extract(i, 'SRS\\d+')
    mat <- read_tsv(i)
	m <- mat[,3, drop = FALSE] %>% as.matrix()
	row.names(m) <- mat$target_id
	colnames(m) <- sample
    mito_sum <- m[mito_genes[mito_genes %in% row.names(m)],1] %>% sum()
	all_sum <- m[,1] %>% sum()
	perc <- (mito_sum/all_sum) * 100
	well_meta[[sample]] <- perc

}



drop_mito <- drop_meta %>% bind_rows() %>% mutate(Barcode = glue::glue("{Barcode}_{sample}")) %>% select(Barcode, `percent.mt`)

well_mito <- well_meta %>% bind_rows() %>% pivot_longer(everything()) %>% arrange(-value)
colnames(well_mito) <- c('Barcode', 'percent.mt')

mito <- bind_rows(well_mito, drop_mito)

write_tsv(mito, file = 'mito_counts.tsv')
