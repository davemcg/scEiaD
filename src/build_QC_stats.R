library(Matrix)
library(Seurat)
library(tidyverse)
library(glue)
args = commandArgs(trailingOnly=TRUE)
git_dir = Sys.getenv('SCIAD_GIT_DIR')

load(args[1])
load(args[2])
#load('umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15.umapFilter.Rdata')

load(args[3])
#load('/data/mcgaugheyd/datashare/scEiaD/2020_08_13/counts_unfiltered.Rdata')

mt <- data.table::fread('pipeline_data/clean_quant/mito_counts.tsv')
mt <- mt %>% unique() %>% group_by(srs, barcode) %>% summarise(`percent.mt` = max(`percent.mt`))
mt <- mt %>% dplyr::rename(value = barcode) %>% mutate(value = gsub(':','_', value))

sobj <- CreateSeuratObject(raw_counts)

#sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
load('pipeline_data/cell_info/cell_info_labelled.Rdata')
meta <- sobj@meta.data %>% as_tibble(rownames = 'value') %>% left_join(cell_info_labels)
meta <- meta %>% left_join(mt )
# add in info for cells removed bc of high mito
meta_study <- meta %>% select(sample_accession, organism, Platform, study_accession) %>% unique()
mt_tossed <- mt[!mt$value %in% meta$value,]  %>% mutate(sample_accession = str_extract(value, '_\\w+') %>% gsub('_','',.))
#meta %>% select(-`percent.mt`) %>% left_join(mt)
mt_tossed <- mt_tossed %>% left_join(., meta_study, by = 'sample_accession') %>% filter(`percent.mt` > 10)



meta <- bind_rows(meta, mt_tossed) %>% filter(!Source %in% c('Cell Culture', 'Organoid'))
# remove excluded samples/studies
exclude <- scan(glue('{git_dir}/data/exclusion.txt'), what = 'character')
cell_info <- data.table::fread('pipeline_data/cell_info/all_cell_info.tsv')
remove_samples <- cell_info %>%  filter(Tissue %in% c('Organoid', 'Cell Culture')) %>% pull(value)
remove_samples <- c(remove_samples, cell_info %>% filter(TissueNote %in% c('Organoid','Cell Culture')) %>% pull(value))
meta <- meta %>% filter(!value %in% remove_samples)
QC <- meta %>% mutate(QC = 
				case_when(value %in% umap$Barcode ~ 'Passed QC',
							nFeature_RNA <= 200 ~ 'Too Few Unique Tx',
							(nFeature_RNA >= 9000) & (Platform %in%
                        c('10xv2','10xv3','DropSeq')) ~ 'Too Many Unique Tx', 
							percent.mt > 10 ~ 'High % MT',
							value %in% umapDoubs$Barcode ~ 'In Silico Doublet',
							TRUE ~ '?'))

write_tsv(QC, path = 'QC.tsv.gz')
