# transfer clark ... blackshaw cell labels to all cells
library(tidyverse)

clark_labels <- read_csv('https://www.dropbox.com/s/y5lho9ifzoktjcs/10x_mouse_retina_development_phenotype.csv?dl=1')
load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/harmony_umap.Rdata')
load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')
sample_table <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv')

samps <- str_extract(harmony_umap %>% row.names(), '(SRS|iPSC_RPE_scRNA_)\\d+')
new <- samps %>% enframe(value = 'sample_accession') %>% left_join(., sra_metadata_extended %>% select(sample_accession, study_accession, Platform, Age, Source, TissueNote) %>% unique())

harmony_umap <- harmony_umap %>% as_tibble(rownames = 'Barcode')
harmony_umap <- cbind(harmony_umap, new)
harmony_umap <- harmony_umap %>% 
  mutate(Age = case_when(Age == 1000 ~ 30, TRUE ~ Age)) %>% 
  left_join(sample_table) %>% 
  as_tibble()
harmony_umap %>% as_tibble()

clark_labels <- clark_labels %>% 
  mutate(core = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>% 
  select(CellType, new_CellType, umap_CellType, umap_coord1, 
         umap_coord2, umap_coord3, X1, barcode, sample, age, 
         core)

harmony_umap_SRP158081 <- harmony_umap %>% mutate(core = gsub('_.*','', Barcode)) %>% 
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

union <- left_join(harmony_umap_SRP158081, clark_labels, 
                   by = c('core', 'sample')) %>% 
  filter(!is.na(CellType))                            

save(union, file = '~/git/massive_integrated_eye_scRNA/data/clark_mcgaughey_label_union.Rdata')
