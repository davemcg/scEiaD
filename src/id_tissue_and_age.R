library(tidyverse)

load('data/sra_metadata.Rdata')

# age is in days, 0 is birth (so you can have negative if there are any prenatal)
# unnest the biosample_attribute_recs for easier parsing later
attribute_info <- list()
for (row in 1:nrow(sra_metadata)){
  #print(row)
  if (!is.null((sra_metadata[row,] %>% pull(biosample_attribute_recs))[[1]])){
    attribute_info[[row]] <- (sra_metadata[row,] %>% 
                                pull(biosample_attribute_recs))[[1]] %>% 
      select(attribute_name, value) %>% 
      spread(attribute_name, value)
  } else {attribute_info[[row]] <- tibble(age = NA)}
}
biosample_attribute_recs <- attribute_info %>% 
  bind_rows() %>% 
  unite(biosample_attribute_recs, remove = FALSE, sep = ' ')
# biosample_attribute_recs <- sra_metadata %>% 
#   select(run_accession, biosample_attribute_recs) %>% 
#   unnest() %>% 
#   group_by(run_accession) %>% 
#   summarise(biosample_attribute_recs = paste(value, collapse = ' ')) 

sra_metadata_extended <- sra_metadata %>% select(-biosample_attribute_recs)
sra_metadata_extended <- cbind(sra_metadata_extended, biosample_attribute_recs)

sra_metadata_extended <- sra_metadata_extended %>% 
  filter(Platform != 'BULK', study_accession != 'SRP149898') %>% 
  mutate(SRP158081_age = str_extract(biosample_title,'^[E|P]\\d+')) %>% 
  mutate(SRP158081_age = case_when(grepl('Nfia', biosample_title) ~ 14,
                                   grepl('^E', SRP158081_age) ~ substr(SRP158081_age, 2,6) %>% as.numeric() - 19,
                                   grepl('^P', SRP158081_age) ~ substr(SRP158081_age, 2,6) %>% as.numeric,
                                   biosample_attribute_recs == 'Retinal Cells P14 Mixed CD1/C57Bl6' ~ 14)) %>% 
  mutate(SRP161678_age = case_when(study_accession == 'SRP161678' ~ str_extract(`time point`, 'Day\\s\\d+') %>% 
                                     gsub('Day ','',.) %>% as.numeric())) %>% 
  mutate(Tissue = 'Retina',
         Source = case_when(study_accession == 'SRP106476' & grepl('hPSC', biosample_attribute_recs) ~ 'hiPSC',
                            study_accession %in% c('SRP136739', 'SRP159286', 'SRP161678') ~ 'hiPSC',
                            study_accession == 'SRP125998' ~ 'Embryonic',
                            study_accession == 'SRP170038' ~ 'Organoid'),  
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
                         study_accession == 'SRP125998' ~ (str_extract(`developmental stage`, '\\d') %>% as.numeric()) * 7,
                         study_accession == 'SRP186396' ~ 19 - (str_extract(`developmental stage`, '\\d+') %>% as.numeric()),
                         study_accession == 'SRP158081' ~ 14,
                         study_accession %in% c('SRP200599','SRP168426') ~ 19 - str_extract(`developmental stage`, '\\d') %>% as.numeric()),   
         Covariate = case_when((study_accession == 'SRP157927' | study_accession == 'SRP158528') & grepl('M1', biosample_title) ~ 'Macaque1',
                               (study_accession == 'SRP157927' | study_accession == 'SRP158528') & grepl('M2', biosample_title) ~ 'Macaque2',
                               (study_accession == 'SRP157927' | study_accession == 'SRP158528') & grepl('M3', biosample_title) ~ 'Macaque3',
                               (study_accession == 'SRP157927' | study_accession == 'SRP158528') & grepl('M4', biosample_title) ~ 'Macaque4',
                               study_accession == 'SRP075719' & grepl('Batch 1', title) ~ 'Batch1',
                               study_accession == 'SRP075719' & grepl('Batch 2', title) ~ 'Batch2',
                               study_accession == 'SRP050054' ~ str_extract(title, 'retina\\s\\d') %>% gsub(' ', '', .),
                               study_accession == 'SRP075720' & grepl('1', `retina id`) ~ 'Batch1',
                               study_accession == 'SRP075720' & grepl('2', `retina id`) ~ 'Batch2',
                               study_accession == 'SRP194595' ~  donor,
                               study_accession == 'SRP075719' ~ sample_accession,
                               study_accession == 'SRP158081' & grepl('rep2', biosample_title) ~ 'Rep2',
                               study_accession == 'SRP158081' & grepl('rep3', biosample_title) ~ 'Rep3',
                               study_accession == 'SRP158081' ~ 'Rep1',
                               study_accession == 'SRP125998' ~ str_extract(title, 'rep\\d'),
                               study_accession == 'SRP168426' ~ str_extract(title, 'E2|F2'),
                               study_accession == 'SRP166660' ~ str_extract(biosample_title, 'run\\d'),
                               study_accession  %in% c('SRP222001','SRP222958') ~ str_extract(title, 'retina \\d') %>% gsub(' ', '', .)),  
         TissueNote = case_when(study_accession == 'SRP158528' | study_accession == 'SRP157927' ~ biosample_title,
                                study_accession == 'SRP073242' ~ 'Vsx2-GFP FACS',
                                study_accession == 'SRP075720' ~ 'Kcng4-cre;stop-YFP X Thy1-stop-YFP Line#1',
                                study_accession == 'SRP106476' ~ 'CRX+/tdTomato OVs and adult retinas were dissociated using papain',
                                study_accession == 'SRP157927' & grepl('M1,', biosample_title) ~ 'Macaque 1, Fovea',
                                study_accession == 'SRP157927' & grepl('M2', biosample_title) ~ 'Macaque 2, Fovea',
                                study_accession == 'SRP157927' & grepl('M3', biosample_title) ~ 'Macaque 3, Fovea',
                                study_accession == 'SRP157927' & grepl('M4', biosample_title) ~ 'Macaque 4, Fovea',
                                study_accession == 'SRP194595' & grepl('Fovea', biosample_title) ~ 'Fovea',
                                study_accession == 'SRP194595' & grepl('Peripheral', biosample_title) ~ 'Peripheral',
                                study_accession == 'SRP158081' ~ 'Sorted Retinal Progenitor Cell (Chx10-GFP positive)',
                                study_accession == 'SRP158528' ~ paste0(source_name, ', ', tissue, ', ', sorting),
                                study_accession == 'SRP159286' ~ 'Retinal organoid',
                                study_accession == 'SRP161678' ~ 'Retinal organoid, human ESC (CRX-GFP line)',
                                study_accession == 'SRP166660' ~ 'Retina, Adult Glast-CreERT2+/tg; Yap5SA+/tg; R26R-tdTomato+/tg mice were induced with tamoxifen followed by NMDA injection. After 2 days, tdTomato+ cells were isolated by fluorescence-activated cell sorting (FACS)',
                                sample_accession == 'SRS4386076' ~ 'Cx3cr1+ single cells from the mouse model of light-induced retinal degeneration',
                                sample_accession == 'SRS4386075' ~ 'Cx3cr1+ single cells, control',
                                sample_accession == 'SRP194595' ~ location,
                                study_accession  %in% c('SRP222001','SRP222958') ~ str_extract(title, 'Macula retina \\d|Peripheral retina \\d'))) %>% 
  dplyr::select(-SRP158081_age, -SRP161678_age)


save(sra_metadata_extended, file = 'data/sra_metadata_extended.Rdata')


