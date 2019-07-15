# run as swarm -f ~/git/massive_integrated_eye_scRNA/src/macaca_10X_bam.swarm --module=aws

# grab the 10X bams for the macaque study
mac_srr <- sra_metadata %>% filter(run_accession %in% (files[(!files %in% existing)] %>% gsub('_.*','',.))) %>% left_join(tech) %>% filter(grepl('Maca', organism)) %>% pull(run_accession)
s3_vector <- vector(mode = 'character', length = length(mac_srr))
for (val in seq(1,length(mac_srr))){
  SRR = mac_srr[val]
  print(SRR)
  page <- GET(paste0('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=', SRR)) %>% 
    content(., 'text', encoding = 'UTF-8')
  split <- (str_split(page, pattern = ' '))[[1]]
  bam_s3 <- tryCatch((split[grepl('bam', split)] %>% 
                        grep("s3",., value = T) %>% 
                        str_split(., '>|<'))[[1]] %>% 
                       grep('s3:',., value = T), error = function(e) return(NA))
  s3_vector[val] <- bam_s3
  #sra_metadata[val, 's3_bam'] <- bam_s3
}


s3_df <- cbind(mac_srr, s3_vector) %>% data.frame()
write(s3_df  %>% rowwise() %>%  mutate(bam = str_split(s3_vector, '/')[[1]][5]) %>% mutate(s3_download = paste0("aws s3 cp ", s3_vector, ' . ; ~/git/massive_integrated_eye_scRNA/src/./bamtofastq ', bam, ' ', mac_srr, ' ; zcat ', mac_srr, '/*/*_R1_001.fastq.gz | gzip > ', mac_srr, '_1.fastq.gz; zcat ', mac_srr, '/*/*_R2_001.fastq.gz | gzip > ', mac_srr, '_2.fastq.gz')) %>% pull(s3_download), file = 'src/macaca_10X_bam.swarm')



######################
# SRP158081 
######################
mac_srr <- sra_metadata %>% filter(run_accession %in% (files[(!files %in% existing)] %>% gsub('_.*','',.))) %>% left_join(tech) %>% filter(study_accession == 'SRP158081') %>% pull(run_accession)
s3_vector <- vector(mode = 'character', length = length(mac_srr))
for (val in seq(1,length(mac_srr))){
  SRR = mac_srr[val]
  print(SRR)
  page <- GET(paste0('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=', SRR)) %>% 
    content(., 'text', encoding = 'UTF-8')
  split <- (str_split(page, pattern = ' '))[[1]]
  bam_s3 <- tryCatch((split[grepl('bam', split)] %>% 
                        grep("s3",., value = T) %>% 
                        str_split(., '>|<'))[[1]] %>% 
                       grep('s3:',., value = T), error = function(e) return(NA))
  s3_vector[val] <- bam_s3
  #sra_metadata[val, 's3_bam'] <- bam_s3
}


s3_df <- cbind(mac_srr, s3_vector) %>% data.frame()
write(s3_df  %>% rowwise() %>%  mutate(bam = str_split(s3_vector, '/')[[1]][5]) %>% mutate(s3_download = paste0("aws s3 cp ", s3_vector, ' . ; ~/git/massive_integrated_eye_scRNA/src/./bamtofastq ', bam, ' ', mac_srr, ' ; zcat ', mac_srr, '/*/*_R1_001.fastq.gz | gzip > ', mac_srr, '_1.fastq.gz; zcat ', mac_srr, '/*/*_R2_001.fastq.gz | gzip > ', mac_srr, '_2.fastq.gz')) %>% pull(s3_download), file = 'src/SRP158081_10X_bam.swarm')






