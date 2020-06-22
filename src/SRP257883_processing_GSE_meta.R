######
# hand downloaded (via safari)
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490473&amp;format=file&amp;file=GSM4490473%5Fdonor%5F22%5Fmac%5FCD31%5Fneg%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490474&amp;format=file&amp;file=GSM4490474%5Fdonor%5F23%5Fmac%5FCD31%5Fpos%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490476&amp;format=file&amp;file=GSM4490476%5Fdonor%5F24%5Fmac%5FCD31%5Fpos%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490479&amp;format=file&amp;file=GSM4490479%5Fdonor%5F25%5Fmac%5FCD31%5Fneg%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490472&amp;format=file&amp;file=GSM4490472%5Fdonor%5F22%5Fmac%5FCD31%5Fpos%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490475&amp;format=file&amp;file=GSM4490475%5Fdonor%5F23%5Fmac%5FCD31%5Fneg%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490478&amp;format=file&amp;file=GSM4490478%5Fdonor%5F25%5Fmac%5FCD31%5Fpos%5Fprocessed%2Etsv%2Egz
# ncbi.nlm.nih.gov/geo/download/?acc=GSM4490477&amp;format=file&amp;file=GSM4490477%5Fdonor%5F24%5Fmac%5FCD31%5Fneg%5Fprocessed%2Etsv%2Egz
##############

library(tidyverse)
files <- list.files('data/', pattern = '.*donor.*mac.*', full.names = TRUE)


file_list <- list()
for (i in files){
  print(i)
  file = basename(i)
  data = read_tsv(i)
  data$sample <- file
  file_list[[i]] <- data %>% select(barcode, donor, celltype, region, age)
}


