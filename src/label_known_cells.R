library(tidyverse)
# load labelled data from clark et all
# https://www.dropbox.com/s/y5lho9ifzoktjcs/10x_mouse_retina_development_phenotype.csv?dl=1
clark_labels <- read_csv('10x_mouse_retina_development_phenotype.csv')

# extract clark blackshaw fields we want
clark_labels <- clark_labels %>% 
  mutate(UMI = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>% 
  select(CellType, new_CellType, umap_CellType, umap_coord1, 
         umap_coord2, umap_coord3, X1, barcode, sample, age, 
         UMI)

## now get macosko labels
# macosko et al
macosko_labels <- read_tsv('http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt', col_names = c('Cell','Cluster'))
macosko_labels <- macosko_labels %>% 
  mutate(CellType = case_when(Cluster == 1 ~ 'Horizontal Cells',
                              Cluster == 2 ~ 'Retinal Ganglion Cells',
                              Cluster < 24 ~ 'Amacrine Cells',
                              Cluster == 24 ~ 'Rods',
                              Cluster == 25 ~ 'Cones',
                              Cluster < 34 ~ 'Bipolar Cells',
                              Cluster == 34 ~ 'Muller Glia',
                              Cluster == 35 ~ 'Astrocytes',
                              Cluster == 36 ~ 'Fibroblasts',
                              Cluster == 37 ~ 'Vascular Endothelium',
                              Cluster == 38 ~ 'Pericytes',
                              Cluster == 39 ~ 'Microglia',
                              TRUE ~ 'None')) %>% 
  mutate(label = gsub('_.*','', Cell),
         UMI = gsub('\\w\\d_','', Cell))

## rgc / bipolar cell labels from karthik shekhar sanes 
# SRP075719
karthik <- read_tsv('~/git/massive_integrated_eye_scRNA/data/shekhar_sanes_ClustAssignFile.txt')
karthik <- karthik %>% 
  # RBC == Rod Bipolar Cell
  mutate(CellType = case_when(CLUSTER == 1 ~ 'Rod Bipolar Cells', 
                              CLUSTER == 2 ~ 'Muller Glia',
                              CLUSTER >= 17 | CLUSTER == 18 | CLUSTER == 19 | CLUSTER == 21 ~ 'Doublet',
                              CLUSTER >= 23 ~ 'Unknown',
                              CLUSTER == 16 ~ 'Amacrine Cells',
                              CLUSTER == 20 ~ 'Rods',
                              CLUSTER == 22 ~ 'Cones', 
                              TRUE ~ 'Bipolar Cells')) %>% 
  mutate(mouse = gsub('_.*', '', CELL_NAME),
         UMI = gsub('.*_','', CELL_NAME)) 
  
## load cell info
load('Mus_musculus_cell_info.Rdata')
# see below for how I got the labelling
cell_info <- cell_info %>% mutate(UMI = gsub('_\\w+', '', value)) %>% 
  mutate(label = case_when(sample_accession == 'SRS866911' ~ 'r2',
                           sample_accession == 'SRS866908' ~ 'r5',
                           sample_accession == 'SRS866912' ~ 'r1',
                           sample_accession == 'SRS866910' ~ 'r3',
                           sample_accession == 'SRS866909' ~ 'r4',
                           sample_accession == 'SRS866907' ~ 'r6',
                           sample_accession == 'SRS866906' ~ 'p1',
                           TRUE ~ 'NOPE'))
meta_SRP158081 <- cell_info %>% 
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
                            TRUE ~ 'UhOh')) %>% 
  left_join(clark_labels %>% 
              select(UMI, sample, new_CellType) %>% 
              dplyr::rename(CellType = new_CellType), 
            by = c('UMI', 'sample')) %>% 
  select(value:batch,CellType) %>% mutate(Paper = 'Clark et al. 2019')

meta_SRP050054 <- cell_info %>% 
  filter(study_accession == 'SRP050054') %>% 
  left_join(macosko_labels %>% 
              select(label, UMI, CellType) , by = c('label', 'UMI' )) %>% 
  select(value:batch, CellType) %>% mutate(Paper = 'Macoko et al. 2015')

meta_SRP075719 <- cell_info %>% 
  filter(study_accession == 'SRP075719') %>% 
  mutate(mouse = case_when(sample_accession == 'SRS1467254' ~ 'Bipolar6',
                           sample_accession == 'SRS1467251' ~ 'Bipolar3',
                           sample_accession == 'SRS1467253' ~ 'Bipolar5',
                           sample_accession == 'SRS1467249' ~ 'Bipolar1',
                           sample_accession == 'SRS1467250' ~ 'Bipolar2',
                           sample_accession == 'SRS1467252' ~ 'Bipolar4',
                           TRUE ~ 'X')) %>% 
  left_join(., karthik %>% select(mouse, UMI, CellType), by = c('UMI', 'mouse')) %>% 
  select(value:batch, CellType) %>% mutate(Paper = 'Shekhar et al. 2016')
  
meta_SRP <- bind_rows(meta_SRP158081, meta_SRP050054, meta_SRP075719)

cell_info_labels <- bind_rows(meta_SRP, 
                              cell_info %>% select(value:batch) %>% 
                                filter(!value %in% meta_SRP$value) %>% 
                                mutate(Paper = NA))

save(cell_info_labels, file = 'Mus_musculus_cell_info_labelled.Rdata')


# ## Figure out what the hell is going on as the above file lists p1, r1 through r6 and GEO lists the 7 samples as 1 - 7....
# ## I'm guessing either goes p1, r1, r2...r6 or r1...r6, p1. 
# ## will use barcode overlap to figure this out
# macosko_barcodes <- colnames(fread('GSE63472_P14Retina_merged_digital_expression.txt.gz', nrows = 1))
# macosko_barcodes <- macosko_barcodes[2:length(macosko_barcodes)]
# ## load "mouse retina 1" and figure out what barcode set it matches best....
# ## SRS866912
# mr1 <- row.names(seurat_late__standard@meta.data) %>% enframe() %>% filter(grepl('SRS866912', value)) %>% mutate(umi = gsub('_.*||\\.\\d*', '', value))
# ## this returns 6450
# ## so retina 1 is r1
# grep('r1_', macosko_barcodes, value = T) %>% gsub('.*_', '', .) %>% enframe() %>% filter(value %in% mr1$umi) %>% dim()
# ## let's see if retina 7 is p1
# ## it is!
# mr7 <- row.names(seurat_late__standard@meta.data) %>% enframe() %>% filter(grepl('SRS866906', value)) %>% mutate(umi = gsub('_.*||\\.\\d*', '', value))
# grep('r7_', macosko_barcodes, value = T) %>% gsub('.*_', '', .) %>% enframe() %>% filter(value %in% mr7$umi) %>% dim()
# ## ok, so appears to p1 == r7
# ## let's just check one in the middle to be sure, r4
# ## yep
# mr4 <- row.names(seurat_late__standard@meta.data) %>% enframe() %>% filter(grepl('SRS866909', value)) %>% mutate(umi = gsub('_.*||\\.\\d*', '', value))
# grep('r4_', macosko_barcodes, value = T) %>% gsub('.*_', '', .) %>% enframe() %>% filter(value %in% mr4$umi) %>% dim()
