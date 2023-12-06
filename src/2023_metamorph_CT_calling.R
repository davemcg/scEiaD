# setup ------
args = commandArgs(trailingOnly=TRUE)
ref <- args[1]
model_type <-  args[2]
cluster_purity <- args[3]

library(tidyverse)
library(metamoRph)
library(uwot)
library(ggplot2)

# load counts
load('/data/mcgaugheyd/datashare/scEiaD/2022_03_22/counts.Rdata')
#load('~/data/scEiaD_2022_02/counts.Rdata')

# load consist diff from plae manuscript
consist_diff <- read_tsv('~/git/eyeMarkers/lists/plae_consist_diff.tsv')
# load metadata for counts
meta_filter <- data.table::fread('/data/mcgaugheyd/projects/nei/mcgaughey/scEiaD_2023_12/metadata_filter.tsv.gz') %>% 
  as_tibble() %>%
  # remove known mis-calls from celltype_predict
  mutate(CellType_predict = case_when(CellType_predict == CellType ~ CellType_predict)) %>% 
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) & !is.na(CellType_predict) ~ 'Tabula Muris',
                                      is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) 

# tweak metadata
meta_filter3 <- meta_filter %>%
  filter(#organism == 'Homo sapiens',
    Organ == 'Eye',
    TechType != 'Well',
    !grepl("RPE|ENDO",sample_accession)) %>%
  mutate(organism = gsub(" ", "_", organism),
         Side = case_when(grepl("Choroid|RPE", Tissue) ~ 'RPEChoroid',
                          grepl("Outf", Tissue) ~ 'OutflowTract',
                          grepl("Iris", Tissue) ~ 'Iris',
                          grepl("Sclera", Tissue) ~ 'Cornea',
                          TRUE ~ Tissue)) %>% 
  mutate(Stage = case_when(as.numeric(Age) <= 10 ~ 'Developing', 
                           as.numeric(Age) > 100 ~ 'Mature',
                           TRUE ~ 'Mature')) %>% 
  mutate(
    ID = paste(Side, Stage, batch, sep = '__'),
    IDq = paste(organism, Side, Stage, study_accession, sep = '__'))

####
# hello darkness ----
####

# cut down to *one* study 
# then roll through each batch within a study and use metamorph/harmony to align with each other
# use the batch with the most labelled cell types (and then cells) as the internal ref
ref <- 'Homo_sapiens__Cornea__Mature__SRP255012'
meta_morph <- meta_filter3 %>%
  filter(IDq == ref)
# pull the batch with the most unique celltypes - break ties on total cell count
batches <- meta_morph %>% 
  group_by(ID, CellType_predict) %>% 
  summarise(Count = n()) %>% 
  summarise(CTCount = n(), Count = sum(Count)) %>% 
  arrange(-CTCount, -Count) %>% pull(ID)
ref_batch <- batches[1]

set.seed(2023-11-15)
# find sd and mean to determine *max* number of cells for a celltype (2 SD)
mm_sd <- meta_morph %>% 
  filter(ID == ref_batch) %>% 
  group_by(CellType_predict) %>% count() %>% 
  pull(n) %>% sd()
mm_mean <- meta_morph %>% 
  filter(ID == ref_batch) %>% 
  group_by(CellType_predict) %>% count() %>%  
  pull(n) %>% mean()
# as a fail-safe for heavy tail dist
# don't let any celltype have more than 1000 
# (many studies have like a bajillion rods)
ref_bcs <- meta_morph %>% 
  filter(ID == ref_batch) %>% 
  group_by(CellType_predict) %>% 
  slice_sample(n = round(min(1000,
                             mm_mean + (2*mm_sd)), 
                         0),
               replace = FALSE) %>% 
  unique() %>% 
  pull(Barcode)
meta_morph_ref <- meta_morph %>% 
  filter(Barcode %in% ref_bcs)

num_pcs <- 20
# run pca with the more balanced cell type representation ----
mm <- run_pca(counts[,meta_morph_ref$Barcode],
              meta_morph_ref,
              ntop = 2000,
              method = "irlba",
              # this abomination add the plae well supported ct diff genes
              # and removes any which have zero total counts
              hvg_force = unique(consist_diff$ID)[!unique(consist_diff$ID) %in% (counts[unique(consist_diff$ID),meta_morph_ref$Barcode] %>% rowSums() %>% enframe() %>% filter(value == 0) %>% pull(name)) ],
              irlba_n = num_pcs)
# loop to use metamorph ----
# to align each batch to the reference
# also uses harmony for moar alignment
morphed_list <- list()
for ( i in batches){
  print(i)
  meta_morph_query <- meta_morph %>% 
    filter(ID == i)
  morphed <- metamoRph(counts[,meta_morph_query$Barcode], 
                       mm$PCA$rotation, mm$center_scale)
  # harmony
  ####
  harmony_meta <- factor(c(rep("new", nrow(morphed)),
                           rep("reference", nrow(mm$PCA$x))),
                         levels = c("reference","new"))
  harmony_embeddings <- harmony::HarmonyMatrix(rbind(morphed, 
                                                     mm$PCA$x),
                                               harmony_meta,
                                               do_pca = FALSE,
                                               reference_values = 'reference',
                                               verbose = FALSE)
  morphed_list[[i]] <- harmony_embeddings[1:nrow(morphed),] %>% data.frame()
}
# new morphed PCs for the study ----
mm_study <- morphed_list %>% bind_rows() 


# do clustering to remove cell type labels ----
### that are in the minority within a cluster
### using a VERY stringent 80% cutoff
g <- bluster::makeSNNGraph(mm_study[,1:num_pcs], k = 20)
clusterN <- igraph::cluster_leiden(g, resolution_parameter = 0.5)$membership
print(unique(clusterN) %>% length())
suspect_bc <- bind_cols(mm_study %>% 
                          as_tibble(rownames = 'Barcode'), tibble(clusterN)) %>% 
  select(Barcode, clusterN) %>% 
  left_join(meta_filter3) %>% group_by(clusterN, CellType_predict) %>% 
  summarise(Count = n(), Barcode = list(Barcode)) %>% 
  mutate(Ratio = Count / sum(Count)) %>% 
  filter(Ratio < cluster_purity) %>% 
  pull(Barcode) %>% 
  unlist()
print(length(suspect_bc))
print(nrow(mm_study))
######################

######################
# build predictor ----
# set Unlabelled CellTypes to NA
ref_bcs <- meta_morph %>%
  filter(!Barcode %in% suspect_bc) %>%  
  group_by(CellType_predict) %>% 
  slice_sample(n = round(min(1000,
                             mm_mean + (2*mm_sd)), 
                         0),
               replace = FALSE) %>% 
  unique() %>% 
  pull(Barcode)
ref_cts <- ref_bcs %>% enframe(value = 'Barcode') %>% 
  left_join(meta_morph) %>% pull(CellType_predict)
ref_cts[ref_cts == 'Unlabelled'] <- NA

model <- model_build(mm_study[ref_bcs,1:num_pcs],
                     ref_cts, model = model_type)
ctp <- model_apply(model,mm_study[,1:num_pcs]) %>% rename(prediction = predict_stringent, Barcode = sample_id)
ctp %>% left_join(meta_filter3) %>% filter(predict == CellType_predict) %>% nrow()
ctp %>% left_join(meta_filter3) %>% filter(predict != CellType_predict) %>% nrow()
#######################

# apply model to all data ---- 
# within the tissue / stage of the reference
ct_info_list <- list()
ct_info_filtered_list <- list()
morphed_list <- list()
for (query in
     meta_morph %>% select(Side, Stage, organism) %>% 
     unique() %>% left_join(meta_filter3, by = c('Side','Stage')) %>% 
     pull(ID) %>% 
     unique()){
  #query <- 'Retina__Mature__SRP222958_DropSeq_retina6'
  print(query)
  meta_morph_query <- meta_filter3 %>% 
    filter(ID == query)
  morphed <- metamoRph(counts[,meta_morph_query$Barcode], mm$PCA$rotation, mm$center_scale)
  harmony_meta <- factor(c(rep("new", nrow(morphed)),
                           rep("reference", nrow(mm$PCA$x))),
                         levels = c("reference","new"))
  harmony_embeddings <- harmony::HarmonyMatrix(rbind(morphed, mm$PCA$x),
                                               harmony_meta,
                                               do_pca = FALSE,
                                               reference_values = 'reference',
                                               max.iter.harmony = 30,
                                               verbose = TRUE)
  morphed <- harmony_embeddings[1:nrow(morphed),] %>% data.frame()
  morphed_list[[query]] <- morphed
  
  
  ct_info <- morphed %>% row.names() %>% enframe(value = 'Barcode') %>% select(Barcode) %>% 
    left_join(meta_filter3 %>% select(Barcode, CellType_predict, CellType))
  ct_info$ref <- ref
  ct_info$predict <-  model_apply(model,morphed[,1:num_pcs]) %>% 
    pull(predict_stringent)
  ct_info$max_score <-  model_apply(model,morphed[,1:num_pcs]) %>% 
    pull(max_score)
  ct_info_list[[query]] <- ct_info
  
  # remove CT level predictions which have terrible accuracy
  # e.g. amacrine cell predictions have <60% accuracy
  filter_ct <- ct_info %>% filter(predict != 'Unknown') %>% 
    group_by(CellType_predict) %>% 
    summarise(accuracy = sum(CellType_predict == predict) / n()) %>% 
    filter(accuracy > 0.6) %>% pull(CellType_predict) %>% unique()
  ct_info_filtered_list[[query]] <- ct_info %>% filter(CellType_predict %in% filter_ct)
}


prediction_info <- ct_info_list %>% 
  bind_rows(.id = 'query') %>% 
  left_join(meta_filter3 %>% select(Barcode, Side, Stage, organism)) %>% 
  relocate(ref,query)

filtered_prediction_info <- ct_info_filtered_list %>% 
  bind_rows(.id = 'query') %>% 
  left_join(meta_filter3 %>% select(Barcode, Side, Stage, organism)) %>% 
  relocate(ref,query)

table(prediction_info %>% filter(Stage == 'Mature', organism == 'Homo_sapiens') %>% pull(CellType_predict), 
      prediction_info %>% filter(Stage == 'Mature', organism == 'Homo_sapiens') %>% pull(predict))

prediction_info %>% filter(Stage == 'Mature', predict == 'Unknown', organism == 'Homo_sapiens') %>% dim()
prediction_info %>% filter(Stage == 'Mature', predict != 'Unknown', organism == 'Homo_sapiens') %>% filter(predict == CellType_predict) %>% dim()
prediction_info %>% filter(Stage == 'Mature', predict != 'Unknown', organism == 'Homo_sapiens') %>% filter(predict != CellType_predict) %>% dim()

morphed_full <- morphed_list %>% bind_rows()

print('one')
write_tsv(prediction_info, file = paste0(gsub(' ','_', ref), "__", model_type, "__", cluster_purity, "__prediction_info.tsv.gz"))
print('two')
write_tsv(filtered_prediction_info, file = paste0(gsub(' ','_', ref), "__", model_type, "__", cluster_purity, "__filtered_prediction_info.tsv.gz"))
save(morphed_full, mm_study, mm, file = paste0(gsub(' ','_', ref), "__", model_type, "__", cluster_purity, ".Rdata"))
