library(tidyverse)

files <- list.files(path = 'pipeline_data/clean_quant', pattern = 'decontX', recursive=TRUE, full.names = TRUE )
files <- files[grep('\\/matrix', files, invert = TRUE)]

spliced_files <- files[grep('_spliced', files)]
unspliced_files <- files[grep('_unspliced', files)]

spliced_l <- list()
for (i in spliced_files){
	spliced_l[[i]] <- read_tsv(i) %>% 
		select(-name) %>% 
		mutate(decontXcontamination = as.numeric(value), 
	    		barcode = as.character(barcode) %>% gsub(':', '_', .)) %>%
		select(-value)
}

unspliced_l <- list()
for (i in unspliced_files){
    unspliced_l[[i]] <- read_tsv(i) %>%
        select(-name) %>%
        mutate(decontXcontamination = as.numeric(value),
                barcode = as.character(barcode) %>% gsub(':', '_', .)) %>%
        select(-value)
}

spliced_decontX_contamination <- spliced_l %>% bind_rows()
unspliced_decontX_contamination <- unspliced_l %>% bind_rows()

write_tsv(spliced_decontX_contamination, file = 'pipeline_data/decontX_contamination/spliced.tsv.gz')
write_tsv(unspliced_decontX_contamination, file = 'pipeline_data/decontX_contamination/unspliced.tsv.gz')
