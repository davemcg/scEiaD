# add biosample to old metadata
library(XML)
library(reutils)
library(tidyverse)

orig_meta <- data.table::fread('~/git/scEiaD/data/sample_run_layout_organism_tech.tsv') %>% as_tibble()

samn_getter <- function(srr){
  (efetch(c(srr), db = 'sra', retmode = 'xml'))$content %>% str_extract(., 'SAMN\\d+')
}

samn_l <- list()
for (i in orig_meta %>% 
     filter(grepl('SRR', run_accession), !run_accession %in% names(samn_l)) %>% 
     pull(run_accession) %>% unique()){
  print(i)
  samn_l[[i]] <- try({samn_getter(i)})
  Sys.sleep(1)
  
}

orig_meta <- orig_meta %>% left_join(., 
                        samn_l %>% as_tibble() %>% 
                          pivot_longer(cols = everything()) %>% 
                          rename(run_accession = name, biosample = value)) %>% 
  mutate(Organ = case_when(Tissue == 'Heart_and_Aorta' ~ 'Heart',
                           Tissue == 'Limb_Muscle' ~ 'Muscle',
                           Tissue %in% c('Organoid','Retina','RPE','RPE-Choroid') ~ 'Eye',
                           TRUE ~ Tissue))

write_tsv(orig_meta, file = '~/git/scEiaD/data/sample_run_layout_organism_tech_biosample_organ_2021_06_05.tsv')
         
         