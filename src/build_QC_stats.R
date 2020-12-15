library(Matrix)
library(Seurat)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

load(args[1])
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umapFilter.Rdata')

load(args[2])
#load('/data/mcgaugheyd/datashare/scEiaD/2020_08_13/counts_unfiltered.Rdata')

mt <- read_tsv('../mito_counts.tsv')
mt <- mt %>% group_by(Barcode) %>% summarise(`percent.mt` = mean(percent.mt))
sobj <- CreateSeuratObject(raw_counts)

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
load('pipeline_data/cell_info/cell_info_labelled.Rdata')
meta <- sobj@meta.data %>% as_tibble(rownames = 'value') %>% left_join(cell_info_labels)
meta <- meta %>% select(-`percent.mt`) %>% left_join(mt %>% rename(value = Barcode))
# add in info for cells removed bc of high mito
meta_study <- meta %>% select(sample_accession, organism, Platform, study_accession) %>% unique()
mt_tossed <- mt[!mt$Barcode %in% meta$value,]  %>% mutate(sample_accession = str_extract(Barcode, '_\\w+') %>% gsub('_','',.))
meta %>% select(-`percent.mt`) %>% left_join(mt %>% rename(value = Barcode))
mt_tossed <- mt_tossed %>% left_join(., meta_study, by = 'sample_accession') %>% filter(`percent.mt` > 10)

meta <- bind_rows(meta, mt_tossed)
QC <- meta %>% mutate(QC = 
				case_when(value %in% umap$Barcode ~ 'Passed QC',
							nFeature_RNA <= 200 ~ 'Too Few Unique Tx',
							(nFeature_RNA >= 3000) & (Platform %in%
                        c('10xv2','10xv3','DropSeq')) ~ 'Too Many Unique Tx', 
							percent.mt > 10 ~ 'High % MT',
							value %in% umapDoubs$Barcode ~ 'In Silico Doublet',
							TRUE ~ 'Passed QC'))

write_tsv(QC, path = 'QC.tsv.gz')
