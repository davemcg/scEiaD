# Rscript

## extract sra project accessions with GEO ID in the dumbest way possible ->
## by building a url and dumping the web page into R

library(httr)
library(tidyverse)
library(DBI)


# function to open page for GSE and extract info
gse_info_maker <- function(gse){  
  # load gse page and get summary and design and prj
  page <- GET(paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', gse)) %>% 
    content(., 'text', encoding = 'UTF-8')
  summary <- ((str_split(page, pattern = '<td nowrap>'))[[1]] %>% 
                grep("Summary<", ., value = TRUE) %>% 
                str_split(., pattern = 'justify\\">|<br>'))[[1]][2]
  design <- ((str_split(page, pattern = '<td nowrap>'))[[1]] %>% 
               grep("Overall design<", ., value = TRUE) %>% 
               str_split(., pattern = 'justify\\">|<br>'))[[1]][2]
  split <- (str_split(page, pattern = ' '))[[1]]
  #acc <- split[grepl('acc', split)]
  acc <- (str_split(split, pattern = '=|>|<|\\"')) %>% unlist()
  prj <- acc[grepl('SRP', acc)] %>% unique() 
  if (length(prj) == 0){prj <- NA}
  line = c(gse, prj, summary, design)
  names(line) <- c('GSE', 'SRA_PROJECT_ID', 'Summary', 'Design')
  line
}
# web search to get info on GEO IDs
base_url <- 'https://www.ncbi.nlm.nih.gov/gds/?term='
GEO_ids <- scan('data/GEO_IDs.txt', what = 'character')
search_url <- paste0(base_url, paste(GEO_ids, collapse = '+'))
gse_prj <- data.frame(matrix(ncol = 4, nrow = 0))
# extract GEO ID and SRA project ID from web page
colnames(gse_prj) <- c("GSE", "SRA_PROJECT_ID", "Summary", "Design")
for (i in GEO_ids){
  search_url <- paste0(base_url, i)
  # load GEO search page to extract GSE
  page <- GET(search_url) %>% httr::content(., type = 'text', encoding = 'UTF-8')
  split <- (str_split(page, pattern = ' '))[[1]]
  acc <- split[grepl('acc', split)]
  acc <- (str_split(acc, pattern = '\\=|\\"|\\[|\\]|\\<|\\>')) %>% unlist()
  

  gse <- acc[grepl('GSE', acc)] %>% unique() 
  line <- gse_info_maker(gse)
  

  gse_prj <- bind_rows(gse_prj, line)
}
# super-series
# some GSE are a bigger collection of other GSE....yay thanks GEO
# if that's the case, then parse out all the related GSE
# and re-run loop
gse_prj_2 <- data.frame(matrix(ncol = 4, nrow = 0))
# extract GEO ID and SRA project ID from web page
colnames(gse_prj_2) <- c("GSE", "SRA_PROJECT_ID", "Summary", "Design")
for (gse in (gse_prj %>% filter(grepl('SuperSeries', Summary)) %>% pull(GSE))){
  page <- GET(paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', gse)) %>% 
    content(., 'text', encoding = 'UTF-8')
  split <- (str_split(page, pattern = ' '))[[1]]
  acc <- split[grepl('acc=GSE', split)]
  acc <- (str_split(acc, pattern = '=|>|<|\\"|&')) %>% unlist()
  gse_sub <- acc[grepl('GSE', acc)] %>% unique() 
  # remove original gse
  gse_sub <- gse_sub[!(gse_sub == gse)]
  for (gse_ID in gse_sub){
    print(gse_ID)
    line <- gse_info_maker(gse_ID)
    print(line)
    gse_prj_2 <- bind_rows(gse_prj_2, line)
  }
}

# missing studies
# hand check to make certain you meant to keep
# on 2019 06 28 the NA were micro-array based studies
# GSE81905 is a super-series and grabbed by gse_prj_2
gse_prj %>% filter(is.na(SRA_PROJECT_ID))
# merge back together
gse_prj <- bind_rows(gse_prj, gse_prj_2) %>% 
  filter(!is.na(SRA_PROJECT_ID), !grepl('bulk RNA', Summary)) %>% # hand remove GSE81902 as this is bulk RNA-seq
  unique()

con <- dbConnect(bigrquery::bigquery(),
                 project = 'isb-cgc-01-0006',
                 dataset = 'omicidx',
                 billing = 'rosy-solstice-244919')

# do not select all columns, as that costs more money
# bigquery essentially charges by data searched
# and if you exclude columns you reduce the amount
# of data you search across

# build query
sql_experiment <- paste0("SELECT study_accession, sample_accession, alias,
                          experiment_accession, title, attributes, instrument_model, 
                          library_layout, library_strategy, library_construction_protocol, 
                          library_layout_length, library_layout_sdev, library_source, platform
              FROM sra_experiment WHERE study_accession IN ('", 
              paste(gse_prj %>% 
                      filter(!is.na(SRA_PROJECT_ID)) %>% 
                               pull(SRA_PROJECT_ID), 
                    collapse = "','"), "')")
# run query
sra_experiment <- dbGetQuery(con, sql_experiment)

sql_sample <- paste0("SELECT sample_accession, organism, taxon_id, BioSample
              FROM sra_sample WHERE accession IN ('", 
                     paste(sra_experiment$sample_accession, 
                           collapse = "','"), "')")
sra_sample <- dbGetQuery(con, sql_sample)

sql_biosample <- paste0("SELECT accession, attributes, attribute_recs, title, taxonomy_name
              FROM biosample WHERE accession IN ('", 
                        paste(sra_sample$BioSample, 
                              collapse = "','"), "')")

sql_run <- paste0("SELECT accession, experiment_accession
              FROM sra_run WHERE experiment_accession IN ('", 
                  paste(sra_experiment$experiment_accession, 
                        collapse = "','"), "')")

sra_biosample <- dbGetQuery(con, sql_biosample)
sra_run <- dbGetQuery(con, sql_run)


# join the 4 into 1
sra_metadata <- sra_experiment %>% 
  left_join(., sra_sample, by = 'sample_accession') %>% 
  left_join(., sra_biosample %>% 
              rename(BioSample = accession, biosample_attributes = attributes, biosample_attribute_recs = attribute_recs, biosample_title = title),
            by = 'BioSample') %>% 
  left_join(., sra_run %>% rename(run_accession = accession),
            by = 'experiment_accession')

# add in tech info about seq
# have to custom correct SRP158081 as they use a blend of smart-seq and 10X....
tech <- read_tsv('data/scTech.tsv')
sra_metadata <- left_join(sra_metadata, tech) %>% 
  mutate(Platform = case_when(study_accession == 'SRP158081' & grepl('Cell', title) ~ 'SMARTSeq_v2',
                              TRUE ~ Platform),
         UMI = case_when(study_accession == 'SRP158081' & grepl('Cell', title) ~ 'NO',
                              TRUE ~ UMI)
         )

save(sra_metadata, file = 'data/sra_metadata.Rdata')
write_tsv(gse_prj %>% rename(study_accession = 'SRA_PROJECT_ID'), path = 'data/GEO_Study_Level_Metadata.tsv')
# hand add our internal RPE samples
core_rpe = data.frame(sample_accession = c(rep('iPSC_RPE_scRNA_01', 24), rep('iPSC_RPE_scRNA_02', 24), rep('iPSC_RPE_scRNA_03', 24)),
             run_accession = c('scRNA_01_S1_L001_HL7H3BCX2','scRNA_01_S1_L001_HLCLYBCX2','scRNA_01_S1_L001_HLFW7BCX2','scRNA_01_S1_L002_HL7H3BCX2','scRNA_01_S1_L002_HLCLYBCX2','scRNA_01_S1_L002_HLFW7BCX2','scRNA_01_S2_L001_HL7H3BCX2','scRNA_01_S2_L001_HLCLYBCX2','scRNA_01_S2_L001_HLFW7BCX2','scRNA_01_S2_L002_HL7H3BCX2','scRNA_01_S2_L002_HLCLYBCX2','scRNA_01_S2_L002_HLFW7BCX2','scRNA_01_S3_L001_HL7H3BCX2','scRNA_01_S3_L001_HLCLYBCX2','scRNA_01_S3_L001_HLFW7BCX2','scRNA_01_S3_L002_HL7H3BCX2','scRNA_01_S3_L002_HLCLYBCX2','scRNA_01_S3_L002_HLFW7BCX2','scRNA_01_S4_L001_HL7H3BCX2','scRNA_01_S4_L001_HLCLYBCX2','scRNA_01_S4_L001_HLFW7BCX2','scRNA_01_S4_L002_HL7H3BCX2','scRNA_01_S4_L002_HLCLYBCX2','scRNA_01_S4_L002_HLFW7BCX2','scRNA_02_S5_L001_HL7H3BCX2','scRNA_02_S5_L001_HLCLYBCX2','scRNA_02_S5_L001_HLFW7BCX2','scRNA_02_S5_L002_HL7H3BCX2','scRNA_02_S5_L002_HLCLYBCX2','scRNA_02_S5_L002_HLFW7BCX2','scRNA_02_S6_L001_HL7H3BCX2','scRNA_02_S6_L001_HLCLYBCX2','scRNA_02_S6_L001_HLFW7BCX2','scRNA_02_S6_L002_HL7H3BCX2','scRNA_02_S6_L002_HLCLYBCX2','scRNA_02_S6_L002_HLFW7BCX2','scRNA_02_S7_L001_HL7H3BCX2','scRNA_02_S7_L001_HLCLYBCX2','scRNA_02_S7_L001_HLFW7BCX2','scRNA_02_S7_L002_HL7H3BCX2','scRNA_02_S7_L002_HLCLYBCX2','scRNA_02_S7_L002_HLFW7BCX2','scRNA_02_S8_L001_HL7H3BCX2','scRNA_02_S8_L001_HLCLYBCX2','scRNA_02_S8_L001_HLFW7BCX2','scRNA_02_S8_L002_HL7H3BCX2','scRNA_02_S8_L002_HLCLYBCX2','scRNA_02_S8_L002_HLFW7BCX2','scRNA_03_S10_L001_HL7H3BCX2','scRNA_03_S10_L001_HLCLYBCX2','scRNA_03_S10_L001_HLFW7BCX2','scRNA_03_S10_L002_HL7H3BCX2','scRNA_03_S10_L002_HLCLYBCX2','scRNA_03_S10_L002_HLFW7BCX2','scRNA_03_S11_L001_HL7H3BCX2','scRNA_03_S11_L001_HLCLYBCX2','scRNA_03_S11_L001_HLFW7BCX2','scRNA_03_S11_L002_HL7H3BCX2','scRNA_03_S11_L002_HLCLYBCX2','scRNA_03_S11_L002_HLFW7BCX2','scRNA_03_S12_L001_HL7H3BCX2','scRNA_03_S12_L001_HLCLYBCX2','scRNA_03_S12_L001_HLFW7BCX2','scRNA_03_S12_L002_HL7H3BCX2','scRNA_03_S12_L002_HLCLYBCX2','scRNA_03_S12_L002_HLFW7BCX2','scRNA_03_S9_L001_HL7H3BCX2','scRNA_03_S9_L001_HLCLYBCX2','scRNA_03_S9_L001_HLFW7BCX2','scRNA_03_S9_L002_HL7H3BCX2','scRNA_03_S9_L002_HLCLYBCX2','scRNA_03_S9_L002_HLFW7BCX2'),
             library_layout = rep('PAIRED', 72),
             organism = rep('Homo sapiens', 72),
             Platform = rep('10xv2', 72),
             UMI = rep('YES', 72),
             study_accession = rep('OGVFB_Hufnagel_iPSC_RPE', 72),
             stringsAsFactors = FALSE)
             
# remove BULK RNA-seq and SRP149898 which is missing the crucial paired end reads (need to contact author)
write_tsv(bind_rows(sra_metadata %>% 
            select(sample_accession, run_accession, library_layout, organism, Platform, UMI, study_accession) %>% 
            filter(Platform != 'BULK', 
                   study_accession != 'SRP149898'), 
            core_rpe),
          path = 'data/sample_run_layout_organism_tech.tsv')
write_tsv(sra_metadata %>% group_by(organism, Platform) %>% sample_n(1) %>% select(sample_accession, run_accession, library_layout, organism, Platform, UMI), path = 'data/sample_run_layout_organism_tech_for_svg.tsv')



      