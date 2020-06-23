############
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037981&amp;format=file&amp;file=GSM4037981%5Fmacula%5Fdonor%5F1%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037982&amp;format=file&amp;file=GSM4037982%5Fmacula%5Fdonor%5F2%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037983&amp;format=file&amp;file=GSM4037983%5Fmacula%5Fdonor%5F3%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037984&amp;format=file&amp;file=GSM4037984%5Fperipheral%5Fdonor%5F1%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037985&amp;format=file&amp;file=GSM4037985%5Fperipheral%5Fdonor%5F2%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037986&amp;format=file&amp;file=GSM4037986%5Fperipheral%5Fdonor%5F3%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037987&amp;format=file&amp;file=GSM4037987%5Fmacula%5Fdonor%5F4%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037988&amp;format=file&amp;file=GSM4037988%5Fmacula%5Fdonor%5F5%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037989&amp;format=file&amp;file=GSM4037989%5Fmacula%5Fdonor%5F6%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037990&amp;format=file&amp;file=GSM4037990%5Fmacula%5Fdonor%5F7%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037991&amp;format=file&amp;file=GSM4037991%5Fperipheral%5Fdonor%5F4%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037992&amp;format=file&amp;file=GSM4037992%5Fperipheral%5Fdonor%5F5%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037993&amp;format=file&amp;file=GSM4037993%5Fperipheral%5Fdonor%5F6%5Fenriched%5Fexpression%2Etsv%2Egz
# www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4037994&amp;format=file&amp;file=GSM4037994%5Fperipheral%5Fdonor%5F7%5Fenriched%5Fexpression%2Etsv%2Egz
###################

library(tidyverse)
files <- list.files('data/', pattern = 'GSM403.*.*.*donor.*tsv.*', full.names = TRUE)


file_list <- list()
for (i in files){
  print(i)
  file = basename(i)
  data = read_delim(i, delim = ' ')
  data$sample <- file
  data <- data %>% select(final_cluster_labels, library, sample, PLP1, NCMAP, SCN7A, MLANA, VWF, ACTA2, IGF2, RPE65, BEST1, PTPRC, CD79A, CD2, AIF1, KIT)
  colnames(data)[1:3] <- c('Barcode', 'Cluster', 'Sample')
  file_list[[i]] <- data %>% select(Barcode, Cluster, Sample, PLP1, NCMAP, SCN7A, MLANA, VWF, ACTA2, IGF2, RPE65, BEST1, PTPRC, CD79A, CD2, AIF1, KIT)
}

SRP218652__non_enriched <- file_list %>% bind_rows() %>% filter(!grepl('enriched', Sample)) %>% 
  mutate(CellType = case_when(Cluster %in% c('1','2') ~ 'Schwann',
                              Cluster == 3 ~ 'Melanocytes',
                              Cluster == 4 ~ 'Endothelial',
                              Cluster == 5 ~ 'Pericytes',
                              Cluster == 6 ~ 'Fibroblasts',
                              Cluster == 7 ~ 'RPE',
                              Cluster == 8 ~ 'B-cell',
                              Cluster == 9 ~ 'T-cell',
                              Cluster == 10 ~ 'Macrophage',
                              Cluster == 11 ~ 'Mast'
                              )) %>% filter(!is.na(CellType))
SRP218652__enriched <- file_list %>% bind_rows() %>% filter(grepl('enriched', Sample)) %>% mutate(CellType = case_when(Cluster %in% c('5','6','7','8') ~ 'Endothelial')) %>%  filter(!is.na(CellType))

SRP218652 <- bind_rows(SRP218652__enriched, SRP218652__non_enriched)
gsm <- str_extract(SRP218652$Sample, 'GSM\\d+') %>% unique()
srs <- c()
for (i in gsm){
  srs <- c(srs, system(paste('/Users/mcgaugheyd/anaconda3/bin/pysradb gsm-to-srs ', i), intern = TRUE)[2] %>% str_extract(., 'SRS\\d+'))
}
conversion <- cbind(gsm, srs) %>% as_tibble()
SRP218652 <- SRP218652 %>% mutate(gsm = str_extract(Sample, 'GSM\\d+')) %>% left_join(., conversion, by = 'gsm') %>% 
  mutate(sample_accession = srs)
save(SRP218652, file = 'data/SRP218652__meta.Rdata')
