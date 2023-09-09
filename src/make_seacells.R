
library(Matrix)
library(tidyverse)

load('/data/mcgaugheyd/projects/nei/mcgaughey/scEiaD_2021_06/pipeline_data/clean_quant/Homo_sapiens/hs-homo_sapiens_full_sparse_matrix.Rdata')

system("mkdir temp_csv; mkdir seacells")

bc <- colnames(all_data)
bc_table <- bc %>% enframe() %>% mutate(sample =  gsub("[ATGC]+_","",value))
# mostly smart-seq  and scrb-seq
n100 <- bc_table %>% group_by(sample) %>% summarise(Count = n()) %>% filter(Count > 99) %>% pull(sample)

for (i in n100){
	print(i)
	csv_file = paste0('temp_csv/', i, '.csv.gz')
    write_csv(all_data[,bc_table %>% filter(sample == i) %>% pull(value)] %>% t() %>% as.matrix()  %>% as_tibble(rownames = 'Barcode'),
                file  = csv_file)
    system(paste("/data/mcgaugheyd/conda/envs/seacells/bin/python ~/git/scEiaD/src/make_seacells.py",
				  csv_file,
				  i,
				  paste0('seacells/', i, '.obs.csv'),
				  paste0('seacells/', i, '.seacell.csv')))
	system(paste0("rm ", csv_file))
}
