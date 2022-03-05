library(tidyverse)
# 2021 05 26 new data
gautan <- read_csv('data/SRP255012_meta.csv')
source('~/git/scEiaD/src/biosample_attribute_grabber.R')
attr <- list()
for (i in gautan$BioSample %>% unique()){
  attr[[i]] <- attribute_df_maker(i)
  Sys.sleep(1)
}

orig_meta <- read_tsv('~/git/scEiaD/data/sample_run_layout_organism_tech_biosample_organ_2022_02_26.tsv')

# [1] "sample_accession"  "run_accession"     "library_layout"    "organism"          "Platform"          "UMI"               "study_accession"   "Tissue"            "Covariate"         "Age"               "integration_group"
# [12] "TissueNote"        "Source"            "bam10x"            "Comment"           "biosample"         "Organ"             "sex"               "retina_region"     "strain"            "biosample_title"   "batch"     

gautanN <- gautan %>% 
  filter(!grepl('RGC', source_name),
         Organism != 'Sus scrofa') %>% 
  mutate(Platform = '10xv3',
         UMI = 'Yes',
         Covariate = NA,
         Age = 'Adult',
         integration_group = 'Late',
         Source = 'Tissue',
         bam10x = NA,
         Comment = Donor,
         Organ = 'Eye',
         sex = NA,
         TissueNote = glue::glue("{tissue} {Donor}")) %>% 
  select(sample_accession = Experiment,
         run_accession = Run,
         library_layout = LibraryLayout,
         organism = Organism,
         Platform,
         UMI,
         study_accession = `SRA Study`,
         Tissue = tissue,
         Covariate,
         Age,
         integration_group,
         TissueNote,
         Source,
         bam10x,
         Comment,
         biosample = BioSample,
         Donor,
         Organ, 
         sex) %>% 
  mutate(Covariate = glue::glue("{Tissue}_{Donor}"),
         Batch = Covariate) %>% 
  left_join(attr %>% bind_rows() %>% filter(attribute == 'biosample_title') %>% rename(biosample = id, biosample_title = value) %>% select(-attribute))



bind_rows(gautanN,orig_meta) %>% head(20) %>% data.frame()
