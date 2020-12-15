args <- commandArgs(trailingOnly = TRUE)
conda_dir =  Sys.getenv('SCIAD_CONDA_DIR')
library(glue)
Sys.setenv(RETICULATE_PYTHON = glue("{conda_dir}/envs/sceasy/bin/python") )
library(sceasy)
library(reticulate)
#loompy <- reticulate::import('loompy')
library(dplyr)

load(args[1]) #umap
load(args[2]) #cluster (meta)
load(args[3]) #seurat_obj
load(args[4]) #cell_info_labels
#load(args[5]) #cell predictions

#colnames(meta) <- c('Barcode', 'cluster')

# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
nmeta <- orig_meta %>% 
  left_join(., cell_info_labels %>% select(-Barcode) %>% select(-contains(c('study_accession', 'Age', 'batch'))) %>% rename(Barcode = value),
			by = 'Barcode') %>%
#  left_join(., predictions %>%
#              as_tibble(rownames = 'Barcode') %>%
#              select(Barcode, CellType_predict = `predicted.id`)) %>%
#  mutate(CellType_predict = case_when(is.na(CellType_predict) ~ CellType,
#                                      TRUE ~ CellType_predict)) %>%
  left_join(meta, by = 'Barcode')
colnames(nmeta)[ncol(nmeta)] <- 'subcluster'
colnames(nmeta)[(ncol(nmeta) - 1)] <- 'cluster'
nmeta <- nmeta %>%  mutate(SubCellType = gsub('p_','', SubCellType)) %>% mutate(SubCellType = gsub('f_','', SubCellType))
nmeta$SubCellType[grepl('^RB|Rods|Peri|^MG$|^Mic$', nmeta$SubCellType)] <- NA

nmeta$CellType[grepl('Astro|Doubl|Fibro|RPE|Vascul|Red ', nmeta$CellType)] <- NA 
nmeta$CellType <- gsub('Early ','',nmeta$CellType)
nmeta$CellType <- gsub('Late  ','',nmeta$CellType)

out <- integrated_obj
out@meta.data <- nmeta %>% data.frame()



sceasy::convertFormat(out, from="seurat", to="anndata",
                       outFile=args[5])
