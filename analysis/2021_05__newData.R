library(tidyverse)
# 2021 05 26 new data
srp292721 <- read_tsv('data/SRP292721_meta.tsv.gz')
srp212788 <- read_tsv('data/SRP212788_meta.tsv.gz')
srp254408 <- read_tsv('data/SRP254408_meta.tsv.gz')
srp255874 <- read_tsv('data/SRP255871_meta.tsv.gz')
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
  select(sample_accession = `EXPERIMENT.STUDY_REF.IDENTIFIERS.PRIMARY_ID`,
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
         Comment) %>% 
  filter(!grepl("_BCR;|_TCR;", TissueNote)) %>% 
  mutate(Organ = case_when(study_accession == 'SRP292721' ~ str_extract(TissueNote, '\\w+_cDNA') %>% gsub('_cDNA','',.),
                           study_accession == 'SRP310237' ~ 'Brain Choroid Plexus',
                           TRUE ~ "Eye"),
         Covariate = case_when(study_accession == 'SRP212788' ~ str_extract(TissueNote, '[\\w\\d]+_') %>% gsub('_R\\d+_','',.)),
         Platform = case_when(study_accession == 'SRP212788' ~ 'SCRBSeq',
                              study_accession == 'SRP292721' ~ '10xv2',
                              study_accession == 'SRP254408' ~ '10xv2'),
         UMI = case_when(study_accession == 'SRP212788' ~ 'No',
                         TRUE ~ 'Yes')
         )


write_tsv(new_meta, '~/git/scEiaD/data/sample_run_layout_organism_tech_2021_newData.tsv')

new_meta %>% 
  filter(sample_accession != )
  mutate(Platform = case_when(study_accession == 'SRP238072') ~ )
