library(tidyverse)

load('data/sra_metadata.Rdata')

# age is in days, 0 is birth (so you can have negative if there are any prenatal)
# unnest the biosample_attribute_recs for easier parsing later
biosample_attribute_recs <- sra_metadata %>% 
  select(run_accession, biosample_attribute_recs) %>% 
  unnest() %>% 
  group_by(run_accession) %>% 
  summarise(biosample_attribute_recs = paste(value, collapse = ' ')) 
sra_metadata_extended <- sra_metadata %>% 
  filter(Platform != 'BULK', study_accession != 'SRP149898') %>% 
  select(-biosample_attribute_recs) %>% 
  left_join(., biosample_attribute_recs) %>% 
  mutate(SRP158081_age = str_extract(biosample_title,'^[E|P]\\d+')) %>% 
  mutate(SRP158081_age = case_when(grepl('^E', SRP158081_age) ~ substr(SRP158081_age, 2,6) %>% as.numeric() - 19,
                                   grepl('^P', SRP158081_age) ~ substr(SRP158081_age, 2,6) %>% as.numeric,
                                   TRUE ~ 1000)) %>% 
  mutate(SRP161678_age = case_when(study_accession == 'SRP161678' ~ str_extract(biosample_attribute_recs, 'Day\\s\\d+') %>% gsub('Day ','',.) %>% as.numeric(),
                                   TRUE ~ 1000)) %>% 
  mutate(Tissue = case_when(study_accession == 'SRP050054' ~ 'Retina',
                            study_accession == 'SRP073242' ~ 'Retina',
                            study_accession == 'SRP075719' ~ 'Retina',
                            study_accession == 'SRP075720' ~ 'Retina',
                            study_accession == 'SRP106476' & !grepl('hPSC', biosample_attribute_recs) ~ 'Retina',
                            study_accession == 'SRP106476' & grepl('hPSC', biosample_attribute_recs) ~ 'Retina',
                            study_accession == 'SRP136739' ~ 'Retina',
                            study_accession == 'SRP157927' ~ 'Retina',
                            study_accession == 'SRP158081' ~ 'Retina',
                            study_accession == 'SRP158528' ~ 'Retina',
                            study_accession == 'SRP159286' ~ 'Retina',
                            study_accession == 'SRP161678' ~ 'Retina',
                            study_accession == 'SRP166660' ~ 'Retina',
                            study_accession == 'SRP186407' ~ 'Retina',
                            study_accession == 'SRP194595' ~ 'Retina'),
         Source = case_when(study_accession == 'SRP106476' & grepl('hPSC', biosample_attribute_recs) ~ 'hiPSC',
                            study_accession %in% c('SRP136739', 'SRP159286', 'SRP161678') ~ 'hiPSC',
                            TRUE ~ 'Tissue'),  
         Age = case_when(study_accession == 'SRP050054' ~ 14,
                         study_accession == 'SRP073242' ~ 17,
                         study_accession == 'SRP075719' ~ 17,
                         study_accession == 'SRP075720' ~ 17,
                         study_accession == 'SRP106476' ~ 90,
                         study_accession == 'SRP157927' ~ 1000,
                         study_accession == 'SRP158081' ~ SRP158081_age,
                         study_accession == 'SRP158528' ~ 1000,
                         study_accession == 'SRP159286' ~ 240,
                         study_accession == 'SRP161678' ~ SRP161678_age,
                         study_accession == 'SRP166660' ~ 1000,
                         study_accession == 'SRP186407' ~ 1000,
                         study_accession == 'SRP194595' ~ 1000,
                         study_accession == 'SRP136739' ~ 90,
                         TRUE ~ 999999),   
         TissueNote = case_when(study_accession == 'SRP073242' ~ 'Vsx2-GFP FACS',
                                study_accession == 'SRP075720' ~ 'Kcng4-cre;stop-YFP X Thy1-stop-YFP Line#1',
                                study_accession == 'SRP106476' ~ 'CRX+/tdTomato OVs and adult retinas were dissociated using papain',
                                study_accession == 'SRP157927' & grepl('M1,', biosample_title) ~ 'Macaque 1, Fovea',
                                study_accession == 'SRP157927' & grepl('M2', biosample_title) ~ 'Macaque 2, Fovea',
                                study_accession == 'SRP157927' & grepl('M3', biosample_title) ~ 'Macaque 3, Fovea',
                                study_accession == 'SRP157927' & grepl('M4', biosample_title) ~ 'Macaque 4, Fovea',
                                study_accession == 'SRP194595' & grepl('Fovea', biosample_title) ~ 'Fovea',
                                study_accession == 'SRP194595' & grepl('Peripheral', biosample_title) ~ 'Peripheral',
                                study_accession == 'SRP158081' ~ 'Sorted Retinal Progenitor Cell (Chx10-GFP positive)',
                                study_accession == 'SRP158528' ~ biosample_attribute_recs %>% gsub('retina ', '', .),
                                study_accession == 'SRP159286' ~ 'Retinal organoid',
                                study_accession == 'SRP161678' ~ 'Retinal organoid, human ESC (CRX-GFP line)',
                                study_accession == 'SRP166660' ~ 'Retina, Adult Glast-CreERT2+/tg; Yap5SA+/tg; R26R-tdTomato+/tg mice were induced with tamoxifen followed by NMDA injection. After 2 days, tdTomato+ cells were isolated by fluorescence-activated cell sorting (FACS)',
                                sample_accession == 'SRS4386076' ~ 'Cx3cr1+ single cells from the mouse model of light-induced retinal degeneration',
                                sample_accession == 'SRS4386075' ~ 'Cx3cr1+ single cells, control',
                                sample_accession == 'SRP194595' ~ str_extract(biosample_attribute_recs, 'Retina.*\\d+') %>% gsub("Retina, ","",.))) %>% 
  dplyr::select(-SRP158081_age, -SRP161678_age)

save(sra_metadata_extended, file = 'data/sra_metadata_extended.Rdata')


