library(tidyverse)
# 2021 05 26 new data
srp292721 <- read_tsv('data/SRP292721_meta.tsv.gz')
srp212788 <- read_tsv('data/SRP212788_meta.tsv.gz')
srp254408 <- read_tsv('data/SRP254408_meta.tsv.gz')
srp255874 <- read_tsv('data/SRP255874_meta.tsv.gz')
srp255871 <- read_tsv('data/SRP255871_meta.tsv.gz')
srp251245 <- read_tsv('data/SRP251245_meta.tsv.gz')
srp238072 <- read_tsv('data/SRP238072_meta.tsv.gz')
srp286543 <- read_tsv('data/SRP286543_meta.tsv.gz')
srp310237 <- read_tsv('data/SRP310237_meta.tsv.gz')

raw_meta <- bind_rows(srp292721,
                      srp212788,
                      srp254408,
                      srp255874,
                      srp255871,
                      srp251245,
                      srp238072,
                      srp286543,
                      srp310237)
orig_meta <- read_tsv('~/git/scEiaD/data/sample_run_layout_organism_tech.tsv')

#  [1] "sample_accession"  "run_accession"     "library_layout"    "organism"         
# [5] "Platform"          "UMI"               "study_accession"   "Tissue"           
# [9] "Covariate"         "Age"               "integration_group" "TissueNote"       
# [13] "Source"            "bam10x"            "Comment"      
new_meta <- raw_meta %>% 
  mutate(Platform = NA,
         UMI = NA,
         Tissue = NA,
         Covariate = NA,
         Age = NA,
         integration_group = NA,
         Source = 'Tissue',
         bam10x = str_extract(bam_url, '\\w+\\.bam'),
         Comment = NA,
         TissueNote = paste0(`EXPERIMENT.TITLE`, ' | ', `SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE`, ' | ', `SAMPLE.TITLE`)) %>% 
  select(sample_accession = `EXPERIMENT.DESIGN.SAMPLE_DESCRIPTOR.IDENTIFIERS.PRIMARY_ID`,
         run_accession = `RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID`,
         library_layout = LibraryLayout,
         organism = ScientificName,
         Platform,
         UMI,
         study_accession = `EXPERIMENT.STUDY_REF.IDENTIFIERS.PRIMARY_ID`,
         Tissue,
         Covariate,
         Age,
         integration_group,
         TissueNote,
         Source,
         bam10x,
         Comment,
         biosample = `SAMPLE.IDENTIFIERS.EXTERNAL_ID`) %>% 
  filter(!grepl("_BCR;|_TCR;", TissueNote),
         run_accession != 'SRR10729793') %>% 
  mutate(Organ = case_when(study_accession == 'SRP292721' ~ str_extract(TissueNote, ':.*[\\w+]+_cDNA') %>% gsub(': |_cDNA.*','',.),
                           study_accession == 'SRP310237' ~ 'Brain',
                           TRUE ~ "Eye"),
         Covariate = case_when(study_accession == 'SRP212788' ~ str_extract(TissueNote, '[\\w\\d]+_') %>% gsub('_R\\d+_','',.),
                               study_accession == 'SRP254408' ~ sample_accession,
                               study_accession == 'SRP255871' ~ str_extract(TissueNote, 'Pt\\d'),
                               study_accession == 'SRP286543' ~ sample_accession,
                               study_accession == 'SRP251245' ~ sample_accession,
                               study_accession == 'SRP310237' ~ paste0(str_extract(bam10x, '^sn|^sc'),
                                                                       '_',
                                                                       str_extract(TissueNote, 'Aged|Adut|Embryo'),
                                                                       '_',
                                                                       str_extract(TissueNote, 'rep\\d+'))),
         Covariate = case_when(study_accession == 'SRP212788' & Covariate == 'M' ~ 'P2',
                               study_accession == 'SRP212788' & Covariate == 'Mac' ~ 'P3',
                               study_accession == 'SRP212788' & Covariate == 'Mac1cell' ~ 'P1',
                               study_accession == 'SRP212788' & Covariate == 'Macu_Nuc' ~ 'P1',
                               study_accession == 'SRP212788' & Covariate == 'Macular' ~ 'P2',
                               study_accession == 'SRP212788' & Covariate == 'Per' ~ 'P3',
                               study_accession == 'SRP212788' & Covariate == 'Periph_Nuc' ~ 'P1' ,
                               study_accession == 'SRP212788' & Covariate == 'sample' ~ 'P3'),
         Platform = case_when(study_accession == 'SRP212788' ~ 'SCRBSeq',
                              study_accession == 'SRP292721' ~ '10xv2',
                              study_accession == 'SRP254408' ~ '10xv2',
                              study_accession == 'SRP255871' ~ '10xv2',
                              study_accession == 'SRP238072' ~ '10xv2',
                              study_accession == 'SRP286543' ~ '10xv2',
                              study_accession == 'SRP255874' ~ '10xv2',
                              study_accession == 'SRP310237' ~ '10xv2',
                              study_accession == 'SRP251245' ~ '10xv3'),
         Tissue = case_when(study_accession == 'SRP251245' ~ 'Outflow Tract',
                            study_accession == 'SRP254408' ~ 'Outflow Tract',
                            study_accession == 'SRP255871' ~ 'Outflow Tract',
                            study_accession == 'SRP255874' ~ 'Outflow Tract',
                            study_accession == 'SRP238072' ~ 'Retina',
                            study_accession == 'SRP310237' ~ 'Brain Choroid Plexus',
                            study_accession == 'SRP238072' ~ 'Retina',
                            study_accession == 'SRP286543' ~ 'Retina',
                            study_accession == 'SRP212788' ~ 'Retina',
                            study_accession == 'SRP292721' ~ Organ,
                            study_accession == 'SRP212788' ~ 'Retina'),
         UMI = case_when(study_accession == 'SRP212788' ~ 'No',
                         TRUE ~ 'Yes'),
         Age = case_when(study_accession == 'SRP292721' ~ 'Adult',
                         study_accession == 'SRP212788' ~ 'Adult',
                         sample_accession == 'SRS6393469' ~ as.character(365 * 73),
                         sample_accession == 'SRS6393471' ~ as.character(365 * 83),
                         sample_accession == 'SRS6393470' ~ as.character(365 * 56),
                         sample_accession == 'SRS6393468' ~ as.character(365 * 55),
                         biosample == 'SAMN14567342' ~ as.character(365 * 65),
                         biosample == 'SAMN14567343' ~ as.character(365 * 56),
                         biosample == 'SAMN14567344' ~ as.character(365 * 64),
                         biosample == 'SAMN14567345' ~ as.character(365 * 60),
                         biosample == 'SAMN14567346' ~ as.character(365 * 78),
                         biosample == 'SAMN14567347' ~ as.character(365 * 78),
                         biosample == 'SAMN14567348' ~ as.character(365 * 74),
                         biosample == 'SAMN14567337' ~ as.character(365 * 10),
                         biosample == 'SAMN14567338' ~ as.character(365 * 10),
                         
                         
         )
  )


write_tsv(new_meta, '~/git/scEiaD/data/sample_run_layout_organism_tech_2021_newData.tsv')

bam_meta <- raw_meta %>% 
  select(`RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID`, bam_url) %>% 
  filter(!is.na(bam_url)) %>% 
  mutate(command = paste0('wget ', bam_url, '; ', 'mv ', str_extract(bam_url, '/\\w+.bam.*') %>% 
                            gsub('/','',.), ' ', `RUN_SET.RUN.IDENTIFIERS.PRIMARY_ID`, '.bam'))

write_tsv(bam_meta, '~/git/scEiaD/data/bam_meta_2021_newData.tsv')
