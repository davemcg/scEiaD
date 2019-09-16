args <- commandArgs(trailingOnly = TRUE)

# transfer clark ... blackshaw cell labels to all cells
library(tidyverse)
library(Seurat)
library(schex)
library(cowplot)
library(future)

plan(strategy = "multicore", workers = 12)
options(future.globals.maxSize = 2400000 * 1024^2)
#load('seurat_merged.Rdata')


# get labelled data from clark et all
clark_labels <- read_csv('https://www.dropbox.com/s/y5lho9ifzoktjcs/10x_mouse_retina_development_phenotype.csv?dl=1')
# macosko et al
macosko_labels <- read_tsv('http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt', col_names = c('Cell','Cluster'))
  # load seurat seurat_merged
load(args[1])
#load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/harmony_umap.Rdata')
load(args[2])
#load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')
sample_table <- read_tsv(args[3])
#sample_table <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv')


meta <- seurat_merged@meta.data %>% as_tibble(rownames = 'Barcode')
# get samples from barcode
samps <- str_extract(meta$Barcode, '(SRS|iPSC_RPE_scRNA_)\\d+')
new <- samps %>% enframe(value = 'sample_accession') %>% 
  # add more info to the metadata
  left_join(., sra_metadata_extended %>% 
              select(sample_accession, study_accession, Platform, Age, Source, TissueNote) %>% 
              unique())


meta <- cbind(meta, new)
meta <- meta %>% 
  mutate(Age = case_when(Age == 1000 ~ 30, TRUE ~ Age)) %>% 
  left_join(sample_table %>% select(-run_accession) %>% unique()) %>% 
  as_tibble()
meta %>% as_tibble()

# extract clark blackshaw fields we want
clark_labels <- clark_labels %>% 
  mutate(core = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>% 
  select(CellType, new_CellType, umap_CellType, umap_coord1, 
         umap_coord2, umap_coord3, X1, barcode, sample, age, 
         core)

# convert our ages and replicate status to match clark blackshaw naming 
# on their sheet
meta_SRP158081 <- meta %>% mutate(core = gsub('_.*','', Barcode)) %>% 
  filter(study_accession == 'SRP158081' & Age != 30 & Platform != 'SMARTSeq_v2') %>% 
  mutate(sample = case_when(Age == -8 ~ 'E11',
                            Age == -7 ~ 'E12_rep1',
                            Age == -5 & Covariate == 'Rep1' ~ 'E14_rep1',
                            Age == -5 & Covariate == 'Rep2' ~ 'E14_rep2',
                            Age == -3 ~ 'E16',
                            Age == -1 & Covariate == 'Rep2' ~ 'E18_rep2',
                            Age == -1 & Covariate == 'Rep3' ~ 'E18_rep3',
                            Age == 0 ~ 'P0',
                            Age == 14 ~ 'P14',
                            Age == 2 & Covariate == 'Rep2' ~ 'P2_rep2',
                            Age == 2 & Covariate == 'Rep3' ~ 'P2_rep3',
                            Age == 5 ~ 'P5',
                            Age == 8 & Covariate == 'Rep1' ~ 'P8_rep1',
                            Age == 8 & Covariate == 'Rep2' ~ 'P8_rep2',
                            TRUE ~ 'UhOh')) 

union <- left_join(meta_SRP158081, clark_labels, 
                   by = c('core', 'sample')) %>% 
  filter(!is.na(CellType))                            

#save(union, file = '~/git/massive_integrated_eye_scRNA/data/clark_mcgaughey_label_union.Rdata')

###########################################################################################
# load('~/git/massive_integrated_eye_scRNA/data/clark_mcgaughey_label_union.Rdata')
# load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')

seurat_merged@meta.data$Age <- meta$Age
seurat_merged@meta.data$Platform <- meta$Platform
seurat_merged@meta.data$study_accession <- meta$study_accession


nmeta <- left_join(meta, union %>% select(Barcode, core:age), by = 'Barcode')
# label the non clark blackshaw labels as missing
nmeta$new_CellType[is.na(nmeta$new_CellType)] <- 'Missing'
seurat_merged@meta.data$new_CellType <- nmeta$new_CellType

DefaultAssay(seurat_merged) <- 'integrated'
missing <- subset(seurat_merged, subset = new_CellType == 'Missing')
labelled <- subset(seurat_merged, subset = new_CellType != 'Missing')
anchors <-  FindTransferAnchors(reference = labelled, 
                                query = missing, 
                                reference.assay = 'integrated',
                                query.assay = 'integrated',
                                dims = 1:30)

predictions <- TransferData(anchorset = anchors, 
                            refdata = labelled@meta.data$new_CellType, 
                            dims = 1:30)

anno <- left_join(nmeta, 
                  predictions %>% as_tibble(rownames = 'Barcode')) %>% 
  mutate(ID = case_when(new_CellType == 'Missing' ~ `predicted.id`, 
                        TRUE ~ new_CellType)) %>% 
  pull(ID) 

meta <- left_join(nmeta, 
                  predictions %>% as_tibble(rownames = 'Barcode')) %>% 
  mutate(ID = case_when(new_CellType == 'Missing' ~ `predicted.id`, 
                        TRUE ~ new_CellType))

seurat_merged@meta.data$new_CellType_transfer <- anno

# save(seurat_merged, file = 'seurat_merged__transfer.Rdata', compress = FALSE)
save(seurat_merged, file = args[4], compress = FALSE)
save(meta, file = args[5])