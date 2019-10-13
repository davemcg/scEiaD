library(tidyverse)
load('data/sra_metadata_extended.Rdata')
# hand add our internal RPE samples
core_rpe = data.frame(sample_accession = c(rep('iPSC_RPE_scRNA_01', 24), rep('iPSC_RPE_scRNA_02', 24), rep('iPSC_RPE_scRNA_03', 24)),
                      run_accession = c('scRNA_01_S1_L001_HL7H3BCX2','scRNA_01_S1_L001_HLCLYBCX2','scRNA_01_S1_L001_HLFW7BCX2','scRNA_01_S1_L002_HL7H3BCX2','scRNA_01_S1_L002_HLCLYBCX2','scRNA_01_S1_L002_HLFW7BCX2','scRNA_01_S2_L001_HL7H3BCX2','scRNA_01_S2_L001_HLCLYBCX2','scRNA_01_S2_L001_HLFW7BCX2','scRNA_01_S2_L002_HL7H3BCX2','scRNA_01_S2_L002_HLCLYBCX2','scRNA_01_S2_L002_HLFW7BCX2','scRNA_01_S3_L001_HL7H3BCX2','scRNA_01_S3_L001_HLCLYBCX2','scRNA_01_S3_L001_HLFW7BCX2','scRNA_01_S3_L002_HL7H3BCX2','scRNA_01_S3_L002_HLCLYBCX2','scRNA_01_S3_L002_HLFW7BCX2','scRNA_01_S4_L001_HL7H3BCX2','scRNA_01_S4_L001_HLCLYBCX2','scRNA_01_S4_L001_HLFW7BCX2','scRNA_01_S4_L002_HL7H3BCX2','scRNA_01_S4_L002_HLCLYBCX2','scRNA_01_S4_L002_HLFW7BCX2','scRNA_02_S5_L001_HL7H3BCX2','scRNA_02_S5_L001_HLCLYBCX2','scRNA_02_S5_L001_HLFW7BCX2','scRNA_02_S5_L002_HL7H3BCX2','scRNA_02_S5_L002_HLCLYBCX2','scRNA_02_S5_L002_HLFW7BCX2','scRNA_02_S6_L001_HL7H3BCX2','scRNA_02_S6_L001_HLCLYBCX2','scRNA_02_S6_L001_HLFW7BCX2','scRNA_02_S6_L002_HL7H3BCX2','scRNA_02_S6_L002_HLCLYBCX2','scRNA_02_S6_L002_HLFW7BCX2','scRNA_02_S7_L001_HL7H3BCX2','scRNA_02_S7_L001_HLCLYBCX2','scRNA_02_S7_L001_HLFW7BCX2','scRNA_02_S7_L002_HL7H3BCX2','scRNA_02_S7_L002_HLCLYBCX2','scRNA_02_S7_L002_HLFW7BCX2','scRNA_02_S8_L001_HL7H3BCX2','scRNA_02_S8_L001_HLCLYBCX2','scRNA_02_S8_L001_HLFW7BCX2','scRNA_02_S8_L002_HL7H3BCX2','scRNA_02_S8_L002_HLCLYBCX2','scRNA_02_S8_L002_HLFW7BCX2','scRNA_03_S10_L001_HL7H3BCX2','scRNA_03_S10_L001_HLCLYBCX2','scRNA_03_S10_L001_HLFW7BCX2','scRNA_03_S10_L002_HL7H3BCX2','scRNA_03_S10_L002_HLCLYBCX2','scRNA_03_S10_L002_HLFW7BCX2','scRNA_03_S11_L001_HL7H3BCX2','scRNA_03_S11_L001_HLCLYBCX2','scRNA_03_S11_L001_HLFW7BCX2','scRNA_03_S11_L002_HL7H3BCX2','scRNA_03_S11_L002_HLCLYBCX2','scRNA_03_S11_L002_HLFW7BCX2','scRNA_03_S12_L001_HL7H3BCX2','scRNA_03_S12_L001_HLCLYBCX2','scRNA_03_S12_L001_HLFW7BCX2','scRNA_03_S12_L002_HL7H3BCX2','scRNA_03_S12_L002_HLCLYBCX2','scRNA_03_S12_L002_HLFW7BCX2','scRNA_03_S9_L001_HL7H3BCX2','scRNA_03_S9_L001_HLCLYBCX2','scRNA_03_S9_L001_HLFW7BCX2','scRNA_03_S9_L002_HL7H3BCX2','scRNA_03_S9_L002_HLCLYBCX2','scRNA_03_S9_L002_HLFW7BCX2'),
                      library_layout = rep('PAIRED', 72),
                      organism = rep('Homo sapiens', 72),
                      Platform = rep('10xv2', 72),
                      UMI = rep('YES', 72),
                      study_accession = rep('OGVFB_Hufnagel_iPSC_RPE', 72),
                      Tissue = c(rep('iPSC', 24), rep('RPE_d42', 24), rep('RPE_Transwell',24)),
                      Covariate = 'None',
                      Age = c(rep(0, 24), rep(42, 24), rep(60, 24)),
                      integration_group = c(rep('iPSC', 24), rep('RPE', 24), rep('RPE',24)),
                      TissueNote = NA,
                      stringsAsFactors = FALSE)

# samples with very low mapped counts
low_n <- c('SRS2728720','SRS2728715','SRS2728718','SRS2728715','SRS2728720','SRS2728718','SRS2728714','SRS2728719','SRS2728714','SRS2728719','SRS2728731','SRS2728694','SRS2728731','SRS2728694','SRS2728685','SRS2728716','SRS2728693','SRS2728685','SRS2728716','SRS2728693','SRS2728730','SRS2728730')

# remove BULK RNA-seq and SRP149898 which is missing the crucial paired end reads (need to contact author)
write_tsv(bind_rows(sra_metadata_extended %>% 
                      select(sample_accession, run_accession, library_layout, organism, Platform, UMI, study_accession, Tissue, Covariate, Age) %>% 
                      mutate(integration_group = case_when(Age < 10  ~ 'Early',
                                               TRUE ~ 'Late')) %>% 
                      filter(Platform != 'BULK', 
                             study_accession != 'SRP149898',
                             !sample_accession %in% low_n) %>% 
                      unique(), 
                    core_rpe),
          path = 'data/sample_run_layout_organism_tech.tsv')
write_tsv(sra_metadata %>% group_by(organism, Platform) %>% sample_n(1) %>% select(sample_accession, run_accession, library_layout, organism, Platform, UMI), path = 'data/sample_run_layout_organism_tech_for_svg.tsv')

