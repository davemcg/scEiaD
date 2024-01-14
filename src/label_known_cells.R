library(tidyverse)
library(yaml)
library(glue)



config=read_yaml(Sys.getenv('SCIAD_CONFIG'))
git_dir=config$git_dir
working_dir=config$working_dir# this is where cell_info lives 

# config <- list()
# git_dir <- '~/git/scEiaD/'
# working_dir= '~/data/scEiaD_2022_02/'
# config$srr_sample_file <- '~/git/scEiaD/data/sample_run_layout_organism_tech_biosample_organ_2022_03_04.tsv'
# config$cell_info <- 'pipeline_data/cell_info/all_cell_info.tsv'
# load(glue('{git_dir}/data/sra_metadata_extended.Rdata'))

setwd(working_dir)
meta <- read_tsv(config$srr_sample_file) %>% select(-TissueNote)

# load labelled data from clark et all
# NO WRONG NOW https://www.dropbox.com/s/y5lho9ifzoktjcs/10x_mouse_retina_development_phenotype.csv?dl=1
clark_labels <- data.table::fread(glue('{git_dir}/data/GSE118614_barcodes.tsv.gz'))
#clark_labels <- read_csv('10x_mouse_retina_development_phenotype.csv')

# extract clark blackshaw fields we1 want
clark_labels <- clark_labels %>% 
  mutate(UMI = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>% 
  select(CellType = umap2_CellType, umap_coord1, 
         umap_coord2, umap_coord3, umap_cluster, V1, barcode, sample, age, 
         UMI)

## now get macosko labels
# macosko et al
# http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt
macosko_labels <- read_tsv(glue('{git_dir}/data/mccaroll_retina_clusteridentities.txt'), col_names = c('Cell','Cluster'))
macosko_labels <- macosko_labels %>% 
  mutate(CellType = case_when(Cluster == 1 ~ 'Horizontal Cells',
                              Cluster == 2 ~ 'Retinal Ganglion Cells',
                              Cluster < 24 ~ 'Amacrine Cells',
                              Cluster == 24 ~ 'Rods',
                              Cluster == 25 ~ 'Cones',
                              Cluster < 34 ~ 'Bipolar Cells',
                              Cluster == 34 ~ 'Muller Glia',
                              Cluster == 35 ~ 'Astrocytes',
                              Cluster == 36 ~ 'Fibroblasts',
                              Cluster == 37 ~ 'Vascular Endothelium',
                              Cluster == 38 ~ 'Pericytes',
                              Cluster == 39 ~ 'Microglia',
                              TRUE ~ 'None')) %>%
  mutate(label = gsub('_.*','', Cell),
         UMI = gsub('\\w\\d_','', Cell))

# wong retina
cell_info <- data.table::fread(config$cell_info) %>% select(-TissueNote)
cell_info <- cell_info %>% filter(study_accession == 'E-MTAB-7316') %>%
				mutate(UMI = gsub('_\\w+', '', value)) %>%
				mutate(sample = case_when(sample_accession == 'ERS2852885' ~ '5',
											sample_accession == 'ERS2852886' ~ '3',
											sample_accession == 'ERS2852887' ~ '4',
											sample_accession == 'ERS2852888' ~ '1',
											sample_accession == 'ERS2852889' ~ '2'))
												
retina_wong <- read_csv(glue('{git_dir}/data/retina_wong_cellbc_cellid.csv')) %>%
				mutate(CellType = gsub(' C\\d+', '', cell.id.orig)) %>% dplyr::select(`cell.bc`, CellType) %>%
				mutate(SubCellType = CellType,
					   CellType = case_when(grepl('RGC', CellType) ~ 'Retinal Ganglion Cells',
											grepl('Amacrine', CellType) ~ 'Amacrine Cells',
											grepl('Bipolar', CellType) ~ 'Bipolar Cells',
											grepl('Cone', CellType) ~ 'Cones',
											grepl('Horizontal', CellType) ~ 'Horizontal Cells',
											grepl('MG', CellType) ~ 'Muller Glia',
											grepl('Rod', CellType) ~ 'Rods',
											TRUE ~ CellType)) %>%
				filter(!is.na(CellType)) %>% 
				filter(CellType != 'Others') %>% 
				separate(`cell.bc`, c('UMI', 'sample'), sep = '-') 

meta_mtab7316 <- cell_info %>% left_join(., retina_wong, by = c('UMI','sample')) %>% mutate(Paper = 'Lukowski et al.')

# sanes mouse rgc crush labels
rgc_crush <- read_tsv(glue('{git_dir}/data/sanes__RGC_Atlas_coordinates.txt'))
rgc_crush <- rgc_crush[-1,]
rgc_crush <- rgc_crush %>% 
				separate(NAME, c('sample', 'UMI'), sep = '_') %>%
				mutate(UMI = gsub('-\\d+','',UMI))
rgc_crush$CellType <- 'Retinal Ganglion Cells'
#rgc_crush$CellType <- NA
## load cell info
cell_info <- data.table::fread(config$cell_info) %>% select(-TissueNote)
cell_info <- cell_info %>% mutate(UMI = gsub('_\\w+', '', value)) %>%
  mutate(sample = case_when(sample_accession == 'SRS5030712' ~ 'aRGC6',
							sample_accession == 'SRS5030713' ~ 'aRGC7',
							sample_accession == 'SRS5030711' ~ 'aRGC5',
							sample_accession == 'SRS5030710' ~ 'aRGC3',
							sample_accession == 'SRS5030714' ~ 'aRGC8',
							sample_accession == 'SRS5030716' ~ 'aRGC10',
							sample_accession == 'SRS5030707' ~ 'aRGC1',
							sample_accession == 'SRS5030708' ~ 'aRGC2',
							sample_accession == 'SRS5030715' ~ 'aRGC9',
							sample_accession == 'SRS5030709' ~ 'aRGC4'))
meta_SRP212151 <- cell_info %>% filter(study_accession == 'SRP212151') %>%
	left_join(., rgc_crush %>% select(sample, UMI, CellType, SubCellType = Cluster), by = c('sample', 'UMI')) %>% 
	select(value:batch,SubCellType, CellType) %>% mutate(Paper = 'Tran et al. 2019')
## rgc / bipolar cell labels from karthik shekhar sanes 
# SRP075719
karthik <- read_tsv(glue('{git_dir}/data/shekhar_sanes_ClustAssignFile.txt'))
karthik <- karthik %>% 
  # RBC == Rod Bipolar Cell
  mutate(CellType = case_when(CLUSTER == 1 ~ 'Rod Bipolar Cells', 
                              CLUSTER == 2 ~ 'Muller Glia',
                              CLUSTER >= 17 | CLUSTER == 18 | CLUSTER == 19 | CLUSTER == 21 ~ 'Doublet',
                              CLUSTER >= 23 ~ 'Unknown',
                              CLUSTER == 16 ~ 'Amacrine Cells',
                              CLUSTER == 20 ~ 'Rods',
                              CLUSTER == 22 ~ 'Cones', 
                              TRUE ~ 'Bipolar Cells')) %>% 
  mutate(SubCellType = case_when(CLUSTER == 7 ~ 'BC1A',
                                 CLUSTER == 9 ~ 'BC1B',
                                 CLUSTER == 10 ~ 'BC2',
                                 CLUSTER == 12 ~ 'BC3A',
                                 CLUSTER == 8 ~ 'BC3B',
                                 CLUSTER == 14 ~ 'BC4',
                                 CLUSTER == 3 ~ 'BC5A',
                                 CLUSTER == 13 ~ 'BC5B',
                                 CLUSTER == 6 ~ 'BC5C',
                                 CLUSTER == 11 ~ 'BC5D',
                                 CLUSTER == 5 ~ 'BC6',
                                 CLUSTER == 4 ~ 'BC7',
                                 CLUSTER == 15 ~ 'BC8/9',
                                 CLUSTER == 1 ~ 'RBC')) %>% 
  mutate(mouse = gsub('_.*', '', CELL_NAME),
         UMI = gsub('.*_','', CELL_NAME)) 

## load cell info
cell_info <- data.table::fread(config$cell_info) %>% select(-TissueNote)
# see below for how I got the labelling
cell_info <- cell_info %>% mutate(UMI = gsub('_\\w+', '', value)) %>% 
  mutate(label = case_when(sample_accession == 'SRS866911' ~ 'r2',
                           sample_accession == 'SRS866908' ~ 'r5',
                           sample_accession == 'SRS866912' ~ 'r1',
                           sample_accession == 'SRS866910' ~ 'r3',
                           sample_accession == 'SRS866909' ~ 'r4',
                           sample_accession == 'SRS866907' ~ 'r6',
                           sample_accession == 'SRS866906' ~ 'p1',
                           TRUE ~ 'NOPE'))
meta_SRP158081 <- cell_info %>% 
  filter(study_accession == 'SRP158081' & Age != 30 & Platform != 'SMARTSeq_v2') %>% 
  mutate(sample = case_when(Age == -8 ~ 'E11',
                            Age == -7 ~ 'E12_rep1',
                            Age == -5 & Covariate == 'Rep1' ~ 'E14_rep1',
                            Age == -5 & Covariate == 'Rep2' ~ 'E14_rep2',
                            Age == -3 ~ 'E16',
                            Age == -1 & Covariate == 'Rep2' ~ 'E18_rep2',
                            Age == -1 & Covariate == 'Rep3' ~ 'E18_rep3',
                            Age == 0 ~ 'P0',
                            Age == 14 ~ 'P14',
                            Age == 2 & Covariate == 'Rep2' ~ 'P2_rep2',
                            Age == 2 & Covariate == 'Rep3' ~ 'P2_rep3',
                            Age == 5 ~ 'P5',
                            Age == 8 & Covariate == 'Rep1' ~ 'P8_rep1',
                            Age == 8 & Covariate == 'Rep2' ~ 'P8_rep2',
                            TRUE ~ 'UhOh')) %>% 
  left_join(clark_labels %>% 
              select(UMI, sample, CellType), 
            by = c('UMI', 'sample')) %>% 
  select(value:batch,CellType) %>% mutate(Paper = 'Clark et al. 2019')



meta_SRP050054 <- cell_info %>% 
  filter(study_accession == 'SRP050054') %>% 
  left_join(macosko_labels %>% 
              select(label, UMI, CellType) , by = c('label', 'UMI' )) %>% 
  select(value:batch, CellType) %>% mutate(Paper = 'Macosko et al. 2015')

meta_SRP075719 <- cell_info %>% 
  filter(study_accession == 'SRP075719') %>% 
  mutate(mouse = case_when(sample_accession == 'SRS1467254' ~ 'Bipolar6',
                           sample_accession == 'SRS1467251' ~ 'Bipolar3',
                           sample_accession == 'SRS1467253' ~ 'Bipolar5',
                           sample_accession == 'SRS1467249' ~ 'Bipolar1',
                           sample_accession == 'SRS1467250' ~ 'Bipolar2',
                           sample_accession == 'SRS1467252' ~ 'Bipolar4',
                           TRUE ~ 'X')) %>% 
  left_join(., karthik %>% select(mouse, UMI, CellType, SubCellType), by = c('UMI', 'mouse')) %>% 
  select(value:batch,CellType, SubCellType) %>% 
  mutate(Paper = 'Shekhar et al. 2016')

# lu clark human dev scRNA
## load cell info
cell_info <- data.table::fread(config$cell_info)  
lu_clark <- data.table::fread(glue('{git_dir}/data/GSE138002_Final_barcodes.csv.gz')) %>%
							mutate(UMI = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>%
				dplyr::rename(CellType = umap2_CellType) %>%
				mutate(SubCellType = CellType, CellType = gsub('BC/Photo_Precurs', 'Photoreceptor Precursors', CellType))
meta_srp223254 <- cell_info %>%
  filter(study_accession %in% c('SRP151023', 'SRP170761','SRP223254')) %>% 
  mutate(UMI = gsub('_\\w+', '', value)) %>% 
  mutate(sample = case_when(Age == -217 ~ '24_Day',
                            Age == -203 ~ 'Hgw11',
                            Age == -196 ~  'Hgw12',
                            Age == -189 ~ 'Hgw13',
                            Age == -182 ~ 'Hgw14',
                            Age == -175 ~ 'Hgw15',
							Age == -168 ~ 'Hgw16',
                            Age == -161 ~ 'Hgw17',
							Age == -154 ~ 'Hgw18',
							sample_accession == 'SRS5434858' ~ 'Hgw19_rep1',
							sample_accession == 'SRS5434858' ~ 'Hgw19_rep2',
							Age == -140 ~ 'Hgw20',
							Age == -126 ~ 'Hgw22',
							sample_accession == 'SRS3443575' ~ 'Hgw24_rep1',
							sample_accession == 'SRS3443577' ~ 'Hgw24_rep2',
							Age == -91 ~ 'Hgw27',
							sample_accession == 'SRS5141025m' ~ 'Hpnd8_rep1',
							sample_accession == 'SRS5141025p' ~ 'Hpnd8_rep2',
							Age == 31360 ~ 'Adult')) %>%
  left_join(lu_clark %>% select(UMI, sample, CellType, SubCellType),
				by = c('UMI', 'sample')) %>% 
  select(value:batch,CellType) %>% mutate(Paper = 'Lu et al. 2020')
## sanes macaque 
cell_info <- data.table::fread(config$cell_info) %>% select(-TissueNote)
sanes_files <- list.files( glue('{git_dir}/data/'), "Macaque*", full.names = TRUE)
sanes <- sanes_files %>% 
  map(read_csv) %>% 
  set_names(sanes_files) %>% 
  bind_rows(.id = 'Type') %>% 
  filter(NAME != 'TYPE') %>% 
  mutate(Type = gsub('_combined.*','', Type)) %>% 
  mutate(Type = gsub('.*_','', Type))

#meta <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv')
meta_MacaqueSanes <- cell_info %>% filter(study_accession %in% c('SRP158528', 'SRP157927')) %>% 
  left_join(sra_metadata_extended %>% 
              select(sample_accession, TissueNote), by = 'sample_accession') %>% 
  mutate(Barcode = gsub("_.*","", value))
meta_MacaqueSanes <- meta_MacaqueSanes %>% 
  left_join(., sanes %>% 
              mutate(NAME = gsub('_S','S', NAME)) %>% 
              separate(NAME, sep = "_|-", 
                       into = c('TissueNote','Barcode','Num'))) %>% 
  mutate(CellType = case_when(Type == 'AC' ~ 'Amacrine Cells',
                              Type == 'BC' ~ 'Bipolar Cells',
                              Type == 'HC' ~ 'Horizontal Cells',
                              Subcluster == 'Endo' ~ 'Vascular Endothelium',
                              Subcluster == 'MG' ~ 'Muller Glia',
                              Subcluster == 'Mic' ~ 'Microglia',
                              Subcluster == 'Pericytes' ~ 'Pericytes',
                              grepl('Cones', Subcluster) ~ 'Cones',
                              Subcluster == 'Rods' ~ 'Rods',
                              Type == 'RGC' ~ 'Retinal Ganglion Cells')) %>% 
  mutate(SubCellType = Cluster) %>% 
  select(value:batch, CellType, SubCellType, TissueNote) %>% 
  mutate(Paper = 'Peng et al. 2019')
# reload
# meta <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv') %>% select(-TissueNote)
# scheetz human fovea
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
scheetz_files <- list.files(glue('{git_dir}/data/'), pattern = 'GSM374599.*donor.*gz', full.names = TRUE)
scheetz <- scheetz_files %>% 
  map(read_tsv, col_types = cols_only(barcode = col_character(), cluster_label = col_character())) %>% 
  set_names(scheetz_files) %>% 
  bind_rows(.id = 'sample') %>% 
  mutate(sample = gsub(".*GSM\\d+_|_expression.*","",sample),
         barcode = gsub('-\\d+','', barcode))
meta_SRP194595 <- cell_info %>% filter(study_accession == 'SRP194595') %>% 
  left_join(., sra_metadata_extended %>% select(sample_accession, biosample_title)) %>% 
  unique() %>% 
  mutate(sample = gsub('\\s+', '_', biosample_title) %>% tolower(),
         barcode = gsub('_.*','', value)) %>% 
  left_join(scheetz) %>% 
  mutate(CellType = case_when(cluster_label %in% c('1','2') ~ 'Rods',
                              cluster_label %in% c('3','4') ~ 'Cones',
                              cluster_label %in% c('5','6') ~ 'Bipolar Cells',
                              cluster_label == '8A' ~ 'Horizontal Cells',
                              cluster_label == '8B' ~ 'Amacrine Cells',
                              cluster_label == '10' ~ 'Pericytes',
                              cluster_label == '11' ~ 'Vascular Endothelium',
                              cluster_label == '12' ~ 'Microglia',
                              cluster_label %in% c('12','13','14','15','16','17') ~ 'Muller Glia')) %>% 
  select(value:batch, CellType) %>% 
  mutate(Paper = 'Voigt et al. Scheetz 2019')

# mennon et al. 
# SRP222001, SRP222958
mennon_seqwell <- read_tsv(glue('{git_dir}/data/GSE137846_Seq-Well_sample_annotations.txt.gz'), 
                           col_names = FALSE, skip = 1) %>% 
  dplyr::select(barcode = X1, Labels = X49)
mennon_drop <- read_tsv(glue('{git_dir}/data/GSE137537_sample_annotations.tsv.gz')) %>% 
  dplyr::select(barcode = Barcode, tissue, Labels)
mennon <- bind_rows(mennon_seqwell, mennon_drop) %>% 
  mutate(SubCelltype = Labels, 
		 CellType = case_when(Labels == 'ACs' ~ 'Amacrine Cells',
                              Labels == 'BPs' ~ 'Bipolar Cells',
                              Labels == 'Endo' ~ 'Vascular Endothelium',
                              Labels == 'HCs' ~ 'Horizontal Cells',
                              Labels == 'Macroglia' ~ 'Muller Glia',
                              Labels == 'RGCs' ~ 'Retinal Ganglion Cells',
                              TRUE ~ Labels),
         Paper = 'Mennon et al. 2019') %>% 
  dplyr::select(-Labels) %>% 
  unique()

meta_mennon <- cell_info %>% filter(study_accession %in% c('SRP222001','SRP222958')) %>% 
  left_join(., sra_metadata_extended %>% select(sample_accession, biosample_title)) %>% 
  unique() %>% 
  mutate(sample = sample_accession,
         barcode = gsub('_.*','', value)) %>% 
  left_join(., mennon %>% 
              separate(barcode, into = c('barcode','sample'), sep = '_|-') %>% 
              mutate(sample = case_when(sample == 'PR8hrs' ~ 'SRS5421628',
                                        sample == 'MR6hrs' ~ 'SRS5421625',
                                        sample == 'MR28hrs' ~ 'SRS5421629',
                                        sample == 'MR8hrs' ~ 'SRS5421627',
                                        sample == 'PR6hrs' ~ 'SRS5421626',
                                        tissue == 'MR' ~ 'SRS5396944',
                                        tissue == 'MR2' ~ 'SRS5396946',
                                        tissue == 'MR3' ~ 'SRS5396948',
                                        tissue == 'PR' ~ 'SRS5396945',
                                        tissue == 'PR2' ~ 'SRS5396947',
                                        tissue == 'PR3' ~ 'SRS5396949')) %>% 
              unique() %>% 
              group_by(barcode, sample) %>% 
              summarise(CellType = paste(CellType, collapse = ',')) %>% 
              filter(!grepl(',', CellType))) %>% 
  mutate(Paper = 'Mennon et al. 2019')

# SRP257883
# choroid/vasculature scheetz/mullins
load(glue('{git_dir}/data/SRP257883_GSE_meta.Rdata'))
SRP257883_GSE_meta <- SRP257883_GSE_meta %>% 
  mutate(SubCellType = celltype, 
         CellType = case_when(celltype == 'smc' ~ 'Smooth Muscle Cell',
                              celltype == 'RBC' ~ 'Red Blood Cells',
                              celltype == 'fibroblast' ~ 'Fibroblasts',
                              celltype == 'macrophage' ~ 'Macrophages',
                              TRUE ~ stringr::str_to_title(celltype))) %>% 
  filter(!CellType %in% c('Rod_rpe')) %>% 
  select(-celltype)
cell_info <- data.table::fread(config$cell_info) %>%  
	filter(study_accession == 'SRP257883') %>%
	mutate(UMI = gsub('_\\w+', '', value)) %>%
	mutate(donor = case_when(sample_accession == 'SRS6517605' ~ 'donor_22_mac_CD31_pos',
								sample_accession == 'SRS6517606' ~ 'donor_22_mac_CD31_neg',
								sample_accession == 'SRS6517607' ~ 'donor_23_mac_CD31_pos',
								sample_accession == 'SRS6517608' ~ 'donor_23_mac_CD31_neg',
								sample_accession == 'SRS6517609' ~ 'donor_24_mac_CD31_pos',
								sample_accession == 'SRS6517610' ~ 'donor_24_mac_CD31_neg',
								sample_accession == 'SRS6517611' ~ 'donor_25_mac_CD31_pos',
								sample_accession == 'SRS6517612' ~ 'donor_25_mac_CD31_neg'))
							
meta_SRP257883 <- cell_info %>% left_join(., SRP257883_GSE_meta %>% 
												mutate(UMI = gsub('_.*','',barcode)), 
											by = c('UMI', 'donor')) %>%
								unique() %>%
								mutate(Paper = 'Voigt et al. 2020')

# tabula muris
cell_info <- data.table::fread(config$cell_info) %>% select(-TissueNote) %>% filter(study_accession == 'SRP131661') 
tm <- data.table::fread(glue('{git_dir}/data/tabula_muris_meta.tsv.gz')) %>%
			dplyr::rename(TabulaMurisCellType = CellType) %>% 
  			mutate(CellType = case_when(TabulaMurisCellType == 'B cell' ~ 'B-Cell',
                              TabulaMurisCellType == 'blood cell' ~ 'Red Blood Cells',
                              TabulaMurisCellType == 'T cell' ~ 'T-Cell',
                              TabulaMurisCellType == 'fibroblast' ~ 'Fibroblasts',
                              grepl('endoth', TabulaMurisCellType, ignore.case = TRUE) ~ 'Endothelial')) %>% 
			select(value, Paper, CellType, TabulaMurisCellType)
meta_TM <- cell_info %>% left_join(., tm, by = 'value' )

## Yan ... Sanes retina rod/cone 
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession == 'SRP255195')
sample_meta <- read_tsv(config$srr_sample_file) %>% filter(study_accession == 'SRP255195')
yan_pr <- read_csv(glue('{git_dir}/data/yan_sanes_PR____Human_retina_combined_all_meta.csv'))
yan_pr <- yan_pr %>% 
  mutate(CellType = case_when(Cluster == 'Astrocytes' ~ 'Astrocytes',
                              grepl('^DB|FMB|IMB|RB1|OFFx|BB', Cluster) ~ 'Bipolar Cells',
                              Cluster == 'Endothelium' ~ 'Endothelial',
                              grepl('Gaba|Gly', Cluster) ~ 'Amacrine Cells',
                              grepl('H1|H2', Cluster) ~ 'Horizontal Cells',
                              grepl('Muller', Cluster) ~ 'Muller Glia',
                              grepl('MG_OFF|MG_ON|RGC|PG_', Cluster) ~ 'Retinal Ganglion Cells',
                              grepl('Cone', Cluster) ~ 'Cones',
                              grepl('Rod', Cluster) ~ 'Rods',
                              Cluster == 'MicroGlia' ~ 'Microglia'),
         SubCellType = Cluster,
         sample = str_extract(NAME, 'H\\d+\\w+S\\d')) %>% 
  filter(!is.na(CellType)) %>% left_join(sample_meta %>%
                                           mutate(sample = str_extract(bam10x, 'H\\d+\\w+S\\d')), 
                                         by = 'sample') %>% 
  filter(!is.na(sample_accession)) %>% 
  separate(NAME, into = c('samp','UMI','lane'), sep = "[_|-]") %>% 
  mutate(value = paste0(UMI, '_', sample_accession)) %>%
  mutate(Paper = 'Yan et al. 2020')
meta_SRP255195 <- cell_info %>% left_join(., yan_pr %>% select(value, Paper, CellType, SubCellType), by = 'value')

# cowan roska
cowan <- data.table::fread(glue('{git_dir}/data/EGAD00001006350_meta.tsv.gz'))
cell_info <- data.table::fread(config$cell_info) %>%
				filter(study_accession == 'EGAD00001006350')

cowan_meta <- cowan %>% 
  mutate(CellType = case_when(grepl('^AC', cell_type)  ~ 'Amacrine Cells',
							  cell_type == 'Ast' ~ 'Astrocytes',
							  cell_type == 'BCell' ~ 'B-Cell',
							  grepl('BC_', cell_type) ~ 'Bipolar Cells',
							  cell_type == 'CM' ~ 'Melanocytes',
							  grepl('END', cell_type) ~ 'Vascular Endothelium',
							  grepl('FB_', cell_type) ~ 'Fibroblasts',
						      grepl('GC_', cell_type) ~ 'Retinal Ganglion Cells',
							  grepl('HC_', cell_type) ~ 'Horizontal Cells',
							  grepl('cone', cell_type) ~ 'Cones',
							  grepl('MAST', cell_type) ~ 'Mast',
						      grepl('MC_', cell_type) ~ 'Muller Glia',
							  grepl('MO_', cell_type) ~ 'Monocyte',
        					  cell_type == 'NK' ~ 'Natural Killer',
							  cell_type == 'PER' ~ 'Pericytes',
							  cell_type == 'RBC' ~ 'Rod Bipolar Cells',
							  cell_type == 'rod' ~ 'Rods',
							  cell_type == 'RPE' ~ 'RPE',
							  cell_type == 'TCell' ~ 'T-Cell',
							  cell_type == 'uG' ~ 'Microglia'), 
		SubCellType = cell_type,
		UMI = str_extract(cell_id, '[ACTG]+-') %>% gsub('-','', .)) %>% 
  filter(!is.na(CellType)) %>%
  mutate(value = paste0(UMI, '_', EGAF)) %>%
  mutate(Paper = 'Cowan et al. 2020')
meta_EGAD00001006350 <- cell_info %>% left_join(., cowan_meta %>% select(value, Paper, CellType, SubCellType), by = 'value')


## voigt mullins rpe
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession == 'SRP218652')
mullins_files <- list.files(glue('{git_dir}/data/'), pattern = 'GSM40379.*gz', full.names = TRUE)
mullins <- mullins_files %>%
  map(read_delim, delim = ' ', col_types = cols_only(barcode = col_character(), final_cluster_labels = col_character(), library = col_character())) %>%
  set_names(mullins_files) %>%
  bind_rows(.id = 'sample') %>%
 # mutate(barcode = final_cluster_labels,
		 mutate(final_cluster_label = final_cluster_labels) %>%
  select(sample, barcode, final_cluster_label, library) %>% 
  mutate(sample = gsub(".*GSM\\d+_|_expression.*","",sample),
         barcode = gsub('-\\d+|_\\d+','', barcode)) %>%
  mutate(sample = paste0(str_extract(sample, 'macula|peripheral'), '_donor', str_extract(sample, '\\d+')))
sample_meta <- read_tsv(config$srr_sample_file) %>% filter(study_accession == 'SRP218652') %>% select(-run_accession) %>% unique()
meta_SRP218652 <- cell_info %>% filter(study_accession == 'SRP218652') %>%
  unique() %>%
  #left_join(sample_meta, by = "sample_accession") %>% 
  mutate(TissueNote = biosample_title) %>% 
  mutate(sample = paste0(retina_region %>% tolower() , '_', 
                         str_extract(batch, 'donor\\d+')),
         barcode = gsub('_.*','', value)) %>%
  left_join(mullins, by =c('sample', 'barcode')) %>%
  mutate(CellType = case_when(grepl('enriched', TissueNote) & final_cluster_label %in% c('1','2') ~ 'Schwann',
                              grepl('enriched', TissueNote) & final_cluster_label %in% c('3') ~ 'Melanocytes',
                              grepl('enriched', TissueNote) & final_cluster_label %in% c('4') ~ 'Endothelial',
                              grepl('enriched', TissueNote) & final_cluster_label == '5' ~ 'Pericytes',
                              grepl('enriched', TissueNote) & final_cluster_label == '6' ~ 'Fibroblasts',
                              grepl('enriched', TissueNote) & final_cluster_label == '8' ~ 'B-Cell',
                              grepl('enriched', TissueNote) & final_cluster_label == '9' ~ 'T-Cell',
							  grepl('enriched', TissueNote) & final_cluster_label == '10' ~ 'Macrophage',
							  grepl('enriched', TissueNote) & final_cluster_label == '11' ~ 'Mast',
                              !grepl('enriched', TissueNote) & final_cluster_label == '1' ~ 'Pericytes',
                              !grepl('enriched', TissueNote) & final_cluster_label == '2' ~ 'Fibroblasts',
                              !grepl('enriched', TissueNote) & final_cluster_label == '3' ~ 'Schwann',
                              !grepl('enriched', TissueNote) & final_cluster_label == '4' ~ 'Fibroblasts',
                              !grepl('enriched', TissueNote) & final_cluster_label %in% c('9','10', '11') ~ 'Macrophage')) %>%
  as_tibble() %>% 
  select(value:batch, CellType) %>%
  unique() %>% 
  mutate(Paper = 'Voigt et al. Mullins 2019')

# SRP259930 Sanes mouse amacrine 63+ cluster paper
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession == 'SRP259930')
metaTN  <- read_tsv(config$srr_sample_file)
sanes_mmAC <- read_csv(glue('{git_dir}/data/MouseAC_metafile.csv'), skip=1) %>% mutate(Mouse = str_extract(TYPE, 'MouseACS\\d+'), Barcode = str_extract(TYPE, '_.*') %>% gsub('_|-\\d+','',.))
meta_SRP259930 <- cell_info %>% 
						mutate(Barcode = gsub("_.*","", value)) %>% 
						select(-TissueNote) %>% 
						#left_join(metaTN %>% select(sample_accession, biosample_title), by = 'sample_accession' ) %>% 
						mutate(Mouse = str_extract(biosample_title, 'MouseACS\\d+')) %>% 
						left_join(sanes_mmAC, by = c('Mouse','Barcode')) %>% 
						mutate(CellType = case_when(!is.na(group) ~ 'Amacrine Cells'), SubCellType = group) %>%
						filter(!is.na(group)) %>%  
						select(value:batch, CellType, SubCellType) %>%
						mutate(Paper = 'Yan et al. 2020')


# heng nathans
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession == 'SRP200499')
heng_files <- list.files(glue('{git_dir}/data/'), pattern = 'GSM38545.*gz', full.names = TRUE)
heng <- heng_files %>%
  map(read_tsv, col_names = FALSE ) %>%
  set_names(heng_files) %>%
  bind_rows(.id = 'sample') %>%
  mutate(Sample = str_extract(sample, "GSM\\d+"),
         Barcode = gsub('-\\d+','', X1),
		 CT = X2) %>%
 select(Sample, Barcode, CT)
meta_SRP200499 <- cell_info %>% 
						mutate(Barcode = gsub("_.*","", value)) %>% 
						#select(-TissueNote) %>% 
						left_join(metaTN %>% select(sample_accession, biosample_title), by = 'sample_accession' ) %>% 
						#mutate(Sample = str_extract(biosample_title, 'GSM38545\\d+')) %>% 
            mutate(Sample = case_when(sample_accession == 'SRX5975001' ~ 'GSM3854512',
                                      sample_accession == 'SRX5975003' ~ 'GSM3854514',
                                      sample_accession == 'SRX5975005' ~ 'GSM3854516',
                                      sample_accession == 'SRX5975007' ~ 'GSM3854518')) %>% 
						left_join(heng, by = c('Sample','Barcode')) %>% 
                        mutate(SubCellType = CT,
							   CellType = case_when(CT == 'Amacrine cells' ~ 'Amacrine Cells',
          							                CT == 'Astrocytes'  ~ 'Astrocytes',
													CT == 'Cone bipolar cells' ~ 'Cone Bipolar Cells',
													CT == 'Rod bipolar cells' ~ 'Rod Bipolar Cells',
													CT == 'Cones' ~ 'Cones',
													CT == 'Horizontal cells' ~ 'Horizontal Cells',
													CT == 'Microglia' ~ 'Microglia',
													CT == 'Monocytes' ~ 'Monocyte',
													CT == 'Muller glia' ~ 'Muller Glia',
													CT == 'Multiplets' ~ 'Doublet',
													CT == 'Retinal ganglion cells' ~ 'Retinal Ganglion Cells',
													CT == 'Rods' ~ 'Rods',
													CT == 'Vascular endothelial cells' ~ 'Endothelial')) %>%
						select(-CT) %>%
						mutate(Paper = 'Heng et al. 2019')




# SRP251245 sanes mouse outflow tract
# SRP255871 sanes human outflow tract
# SRP255874 sanes macaque outflow tract
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession %in%  c('SRP251245', 'SRP255871', 'SRP255874')) %>%
	mutate(sample = case_when(biosample_title == 'MouseTMS2' ~ 'MouseTMS2',
								biosample_title == 'Mulatta2TMR' ~ 'MM1RTM',
								grepl('Pt\\d', biosample_title) ~ str_extract(biosample_title, 'H\\w+'),
								biosample_title == 'Mulatta1TML' ~ 'MM1LTM',
								biosample_title == 'MouseTMS1' ~ 'MouseTMS1'))

outflow <- read_csv(glue("{git_dir}/data/SRP251245_all_five_species_metafile.csv")) %>% filter(NAME != 'TYPE') %>%
				separate(NAME, c('sample','BC'), sep = '_') %>% 
  				mutate(BC = gsub('-\\d+','',BC))  
	
outflow_meta <- cell_info %>% 
					mutate(BC = gsub("_.*","", value)) %>%
					left_join(outflow %>% select(sample, BC, Cluster), by = c('sample','BC')) %>%
					mutate(Cluster = gsub('^\\d+_','',Cluster),
							SubCellType = Cluster) %>% 
					mutate(CellType = case_when(grepl('B Cell|BCell', Cluster, ignore.case = TRUE) ~ 'B-Cell',
													grepl('Beam', Cluster, ignore.case=TRUE) ~ 'Beam',
													grepl('Muscle',Cluster, ignore.case=TRUE) ~ 'Ciliary Muscle',
													grepl('Cornea', Cluster, ignore.case =TRUE) ~ 'Cornea',
													Cluster %in% c('Endothelium', 'JCT','Macrophage',
																	'Melanocyte', 'Pericyte', 'Uveal')  ~ Cluster,
													grepl('JCT', Cluster) ~ 'JCT',
													grepl('Mast', Cluster) ~ 'Mast',
													grepl('schwann', Cluster, ignore.case = TRUE) ~ 'Schwann',
													grepl('NK/T|NKT', Cluster) ~ 'T-Cell',
													grepl('epithe', Cluster) ~ 'Epithelium',
													Cluster == 'SChlemm\'s Canal' ~ 'Schlemm\'s Canal',
													grepl('Endo', Cluster) ~ 'Endothelium')) %>%
					select(-sample, -BC, -Cluster) %>%
					mutate(Paper = 'van Zyl et al. 2020')	
												
#  SRP286543 sanes chick atlas
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession %in%  c('SRP286543'))
chick_atlas <- read_csv(glue('{git_dir}/data/yamagata_sanes_Chick_retina_atlas_meta.csv.gz')) %>% filter(NAME != 'TYPE') %>% 
  separate(NAME, c('Chick','BC'), sep = '_') %>% 
  mutate(BC = gsub('-\\d+','',BC),
         sample_accession = case_when(Chick == 'Chicken1A' ~  'SRS7483320',
							Chick == 'Chicken1B' ~  'SRS7483321',
							Chick == 'Chicken1C' ~  'SRS7483322',
							Chick == 'Chicken1D' ~  'SRS7483323',
							Chick == 'ChickendRGC1' ~  'SRS7483324',
							Chick == 'ChickendRGC2' ~  'SRS7483325',
							Chick == 'ChickenvRGC1' ~  'SRS7483326',
							Chick == 'ChickenvRGC2' ~  'SRS7483327')) 
chick_meta <- chick_atlas %>% select(Chick, BC, Cluster, sample_accession) %>% 
  right_join(cell_info %>% mutate(BC = gsub("_.*","", value)), by = c('BC','sample_accession')) %>% 
  mutate(CellType = case_when(grepl('AC', Cluster) ~ 'Amacrine Cells',
                              grepl('BP', Cluster) ~ 'Bipolar Cells',
                              grepl('HC', Cluster) ~ 'Horizontal Cells',
                              grepl('OG', Cluster) ~ 'Oligodendrocytes',
                              grepl('MG', Cluster) ~ 'Muller Glia',
                              grepl('Cones', Cluster) ~ 'Cones',
                              grepl('Rods', Cluster) ~ 'Rods',
                              grepl('RGC', Cluster) ~ 'Retinal Ganglion Cells'),
		SubCellType = Cluster) %>%
  select(-Chick, -BC, -Cluster) %>%
  mutate(Paper = 'Yamagata et al. 2021')

# SRP292721 human pan atlas
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession %in%  c('SRP292721'))
pan_human <- read_tsv(glue('{git_dir}/data/SRP292721_Annotation_AHCA_alltissues_meta.data_84363_cell.txt.gz')) %>%  
	separate(`...1`, c('Organ', 'cDNA', 'BC'), sep = '_') %>%
	mutate(BC = gsub('-\\d+','',BC)) %>%
    mutate(CellType = case_when(Cell_type_in_merged_data == 'Absorptive Cell' ~ 'Absorptive Cell',
								grepl('B Cell', Cell_type_in_merged_data, ignore.case = TRUE) ~ 'B-Cell',
							 	grepl('Basal Cell', Cell_type_in_merged_data) ~ 'Basal Cell',
								Cell_type_in_merged_data == 'Cholangiocyte' ~ 'Cholangiocyte',
								grepl('Endothelial', Cell_type_in_merged_data) ~ 'Endothelial',
								grepl('Enterocyte', Cell_type_in_merged_data) ~ 'Enterocyte',
								grepl('Fibroblast', Cell_type_in_merged_data) ~ 'Fibroblasts',
								grepl('Erythrocyte', Cell_type_in_merged_data) ~ 'Erthyrocyte',
								grepl('Keratinocyte', Cell_type_in_merged_data) ~ 'Keratinocyte',
								grepl('Macrophage', Cell_type_in_merged_data) ~ 'Macrophage',
								grepl('Melanocyte', Cell_type_in_merged_data) ~ 'Melanocytes',
								grepl('Epithelial', Cell_type_in_merged_data) ~ 'Epithelial',
								grepl('Plasma', Cell_type_in_merged_data) ~ 'Plasma Cell',
								grepl('NK/T|T Cell', Cell_type_in_merged_data) ~ 'T-Cell',
								grepl('Monocyte', Cell_type_in_merged_data) ~ 'Monocyte',
								grepl('Satellite', Cell_type_in_merged_data) ~ 'Satellite Cell',
								grepl('Secretory', Cell_type_in_merged_data) ~ 'Secretory Cell',
								grepl('Smooth Muscle Cell', Cell_type_in_merged_data) ~ 'Smooth Muscle Cell',
								grepl('Tuft Cell', Cell_type_in_merged_data) ~ 'Tuft Cell'))

pan_human_meta <- cell_info %>%
					mutate(BC = gsub("_.*","", value)) %>%
					left_join(pan_human %>% select(Organ, BC, CellType, SubCellType = Cell_type_in_merged_data), by = c('BC','Organ')) %>%
					mutate(Paper = 'He et al. 2020')
					
# SRP310237 mouse brain choroid plexus
cell_info <- data.table::fread(config$cell_info) %>%
    filter(study_accession %in%  c('SRP310237'))
brainCPnuc <- read_tsv(glue('{git_dir}/data/SRP310237__single_nucleus_metadata.txt.gz')) %>% filter(NAME != 'TYPE') %>%
					separate(NAME, c('sample','BC','index'), sep = '\\.') %>%
					mutate(sample_accession = case_when(sample == 'Adult_TCP_1203' ~ 'SRS8433494',
														sample == 'Aged_DCP_1008' ~ 'SRS8433501',
														sample == 'Embryo_HCP_1209' ~ 'SRS8433491',
														sample == 'Aged_HCP_1203' ~ 'SRS8433504',
														sample == 'Adult_HCP_1008' ~ 'SRS8433497',
														sample == 'Embryo_DCP_1008' ~ 'SRS8433493',
														sample == 'Aged_TCP_1203' ~ 'SRS8433499',
														sample == 'Adult_TCP_1008' ~ 'SRS8433492',
														sample == 'Adult_DCP_1203' ~ 'SRS8433496',
														sample == 'Aged_DCP_1203' ~ 'SRS8433502',
														sample == 'Aged_TCP_1008' ~ 'SRS8433498',
														sample == 'Aged_HCP_1008' ~ 'SRS8433505',
														sample == 'Aged_TCP_1209' ~ 'SRS8433500')) %>%
					select(sample_accession, BC, cell_type__ontology_label, cell_type_subset__ontology_label) %>% 
					mutate(CellType = case_when(grepl('Endo', cell_type__ontology_label, ignore.case = TRUE) ~ 'Endothelial',
												grepl('Epith', cell_type__ontology_label, ignore.case = TRUE) ~ 'Epithelial',
												grepl('neuron', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Neuron',
												grepl('Fibroblast', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Fibroblasts',
												grepl('Lymphocyte', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Lymphocyte',
												grepl('^Mac', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Macrophage',
												grepl('Meningeal', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Meningeal',
												grepl('Neutrophil', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Neutrophil',
												grepl('Pericytes', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'Pericytes',
												grepl('SMC', cell_type_subset__ontology_label, ignore.case = TRUE) ~ 'SMC'),
							SubCellType = case_when(!grepl('mesenchymal', cell_type__ontology_label) ~ cell_type__ontology_label,
													TRUE ~ cell_type_subset__ontology_label)) 

	
meta_SRP310237 <- cell_info %>% 
					mutate(BC = gsub("_.*","", value)) %>%
					left_join(brainCPnuc %>% select(sample_accession, BC, CellType, SubCellType), by = c('BC','sample_accession')) %>%
					mutate(Paper = 'Dani et al. 2021')			

# SRP275814 collin lako cornea
# 10pcw
cornea <- data.table::fread(config$cell_info) %>%
    filter(study_accession == 'SRP275814')
pcw10 <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_10PCW_meta.tsv'))
cornea_pcw10 <- cornea %>% filter(grepl('10PCW', biosample_title))
SRP275814_cornea_pcw10 <- cornea_pcw10 %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>% left_join(pcw10) 

cornea_pcw12 <- cornea %>% filter(grepl('12PCW', biosample_title))
pcw12 <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_12PCW_meta.tsv'))
SRP275814_cornea_pcw12 <- cornea_pcw12 %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>% 
							left_join(pcw12 %>% mutate(run = str_extract(cellName, '\\d+'), 
														cellName = str_extract(cellName, '^[ACGTN]+'),
														sample_accession = case_when(run == '2' ~ 'SRX8884955',
																					 run == '1' ~ 'SRX8884954',
																					 run == '4' ~ 'SRX8884957',
																					 run == '3' ~ 'SRX8884956')))
												
cornea_pcw1314 <- cornea %>% filter(grepl('13PCW|14PCW', biosample_title))
pcw1314 <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_13_14PCW_meta.tsv'))
SRP275814_cornea_pcw1314 <- cornea_pcw1314 %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>%
                            left_join(pcw1314 %>% mutate(run = str_extract(cellName, '\\d+'),
                                                        cellName = str_extract(cellName, '^[ACGTN]+'),
                                                        sample_accession = case_when(run == '1' ~ 'SRX8884958',
                                                                                     run == '4' ~ 'SRX8884961',
                                                                                     run == '3' ~ 'SRX8884960',
                                                                                     run == '2' ~ 'SRX8884959')))

cornea_pcw16 <-  cornea %>% filter(grepl('16PCW', biosample_title))
pcw16 <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_16PCW_meta.tsv'))
SRP275814_cornea_pcw16 <- cornea_pcw16 %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>%
                            left_join(pcw16 %>% mutate(run = str_extract(cellName, '\\d+'),
                                                        cellName = str_extract(cellName, '^[ACGTN]+'),
                                                        sample_accession = case_when(run == '1' ~ 'SRX8884962',
                                                                                     run == '2' ~ 'SRX8884963')))

cornea_pcw1718 <- cornea %>% filter(grepl('17PCW|18PCW', biosample_title))
pcw1718 <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_17_18PCW_meta.tsv'))
SRP275814_cornea_pcw1718 <- cornea_pcw1718 %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>%
                            left_join(pcw1718 %>% mutate(run = str_extract(cellName, '\\d+'),
                                                        cellName = str_extract(cellName, '^[ACGTN]+'),
                                                        sample_accession = case_when(run == '1' ~ 'SRX8884964',
                                                                                     run == '2' ~ 'SRX8884965')))
cornea_pcw2021 <- cornea %>% filter(grepl('20PCW|21PCW', biosample_title))
pcw2021 <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_20_21PCW_meta.tsv'))
SRP275814_cornea_pcw2021 <- cornea_pcw2021 %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>%
                            left_join(pcw2021 %>% mutate(run = str_extract(cellName, '\\d+'),
                                                        cellName = str_extract(cellName, '^[ACGTN]+'),
                                                        sample_accession = case_when(run == '1' ~ 'SRX8884966',
                                                                                     run == '2' ~ 'SRX8884967',
																					 run =='3' ~ 'SRX8884968')))
cornea_adult <- cornea %>% filter(!grepl('PCW', biosample_title))
adult <- read_tsv(glue::glue('{git_dir}/data/lako_cornea_adult_meta.tsv'))
SRP275814_cornea_adult <- cornea_adult %>% mutate(cellName = str_extract(value, '^[ANGTC]+')) %>%
                            left_join(adult %>% mutate(run = str_extract(cellName, '\\d+'),
                                                        cellName = str_extract(cellName, '^[ACGTN]+'),
                                                        sample_accession = case_when(run == '1' ~ 'SRX8884969',
                                                                                     run == '2' ~ 'SRX8884971',
																					 run == '4' ~ 'SRX8884972')))
SRP275814 <- bind_rows(SRP275814_cornea_pcw10,
						SRP275814_cornea_pcw12,
						SRP275814_cornea_pcw1314,
						SRP275814_cornea_pcw16,
						SRP275814_cornea_pcw1718,
						SRP275814_cornea_pcw2021,
						SRP275814_cornea_adult) %>% 
							mutate(CellType = case_when(cluster == 'neural crest' ~ 'Neural Crest',
														grepl('corneal epithelium', cluster) ~ 'Corneal Epithelium',
														grepl('corneal endothelium', cluster) ~ 'Corneal Endothelium',
														cluster == 'conjunctival epithelium' ~ 'Conjunctival Epithelium',
														grepl('proliferating', cluster) ~ 'Proliferating Cornea',
														cluster == 'corneal nerves' ~ 'Corneal Nerve',
														grepl('keratocytes', cluster) ~ 'Keratocyte',
														cluster == 'ciliary margin' ~ 'Ciliary Margin',
														grepl('vessel', cluster) ~ 'Blood Vessel',
														grepl('red blood', cluster) ~ 'Red Blood Cell',
														grepl('fibroblast', cluster) ~ 'Fibroblast',
														grepl('mesoderm', cluster) ~ 'Mesoderm',
														cluster == 'melanocytes' ~ 'Melanocyte',
														cluster == 'epithelial basement membrane' ~ 'Corneal Basement Membrane',
														grepl('progenitors', cluster, ignore.case = TRUE) ~ 'Corneal Progenitor',
														cluster == 'limbal stem/progenitor cells' ~ 'Limbal Progenitor', 
														grepl('limbal', cluster, ignore.case = TRUE) ~ 'Limbal'	),
									SubCellType = cluster) %>%
							select(-cluster, -cellName) %>% 
							mutate(Paper = 'Collin et al. 2021')


# Bala
cell_info <- data.table::fread(config$cell_info) %>% select(-TissueNote)
SRP228556 <- read_tsv(glue('{git_dir}/data/Balasubramanian_Zhang_metaData6.tsv.gz')) %>%
				mutate(Paper = 'Balasubramanian et al. 2022', sample_accession = 'SRX7099188') %>% 
				mutate(NAME = gsub('-1','',NAME),
						Barcode = glue('{NAME}_{sample_accession}'),
						Comment = condition,
						SubCellType = CellType,
						CellType = case_when(SubCellType == 'Photoreceptors' ~ 'Photoreceptor Precursor',
											 SubCellType == 'ACs/HCs' ~  'AC/HC_Precursor',
											 SubCellType ==  'RGCs' ~  'Retinal Ganglion Cell',
											 grepl('RPC', SubCellType) ~ 'RPC',
											 grepl('CM', SubCellType) ~ 'Ciliary Margin',
											 grepl('Neurogenic', SubCellType) ~ 'Neurogenic Cell')) %>%
						select(Barcode, CellType, SubCellType, Paper, Comment)

SRP228556_meta <- cell_info %>% filter(sample_accession == 'SRX7099188') %>%
					left_join(SRP228556, by = c('value' = 'Barcode'))
											 
# Gautam Loh
swapper <- cell_info  %>% filter(study_accession == 'SRP255012') %>% select(Covariate, sample_accession)  %>% mutate(samp = gsub('_','',Covariate)) %>% unique()
SRP255012 <- data.table::fread(glue('{git_dir}/data/Gautam_Loh_scp_metadata.txt.gz')) %>% as_tibble() %>% 
					mutate(SubCellType = Sub_celltype,
							CellType = case_when(grepl('T cells', SubCellType) ~ 'T/NK-Cell',
										SubCellType == 'Amacrine cells' ~ 'Amacrine Cell',
										grepl('endothelial', SubCellType) ~ 'Endothelial',
										grepl('ciliary body', SubCellType, ignore.case = TRUE) ~ 'Ciliary Body',
										SubCellType == 'Cone bipolar cells' ~ 'Bipolar Cell',
										SubCellType == 'Cones' ~ 'Cone',
										SubCellType == 'Conjunctival cells' ~ 'Conjunctival Epithelial',
										SubCellType == 'Corneal epithelial cells' ~ 'Corneal Epithelial',
										grepl('fibroblast', SubCellType, ignore.case = TRUE) ~ 'Fibroblast',
										SubCellType == 'Endothelial cells' ~ 'Endothelial',
										SubCellType == 'Horizontal cells' ~ 'Horizontal Cell',
										SubCellType == 'Melanocytes' ~ 'Melanocyte',
										SubCellType == 'Microglia' ~ 'Microglia',
										SubCellType == 'Monocytes' ~ 'Monocyte',
										SubCellType == 'Muller glia' ~ 'Muller Glia',
										grepl('bipolar c', SubCellType) ~ 'Bipolar Cell',
										SubCellType == 'RGCs' ~ 'Retinal Ganglion Cell',
										SubCellType == 'Rod bipolar cells' ~ 'Rod Bipolar Cell',
										SubCellType == 'Rods' ~ 'Rod',
										SubCellType == 'RPE' ~ 'RPE',
										SubCellType == 'Schwann cells' ~ 'Schwann',
										SubCellType == 'Smooth muscle cells' ~ 'Smooth Muscle Cell',
						)) %>%
					left_join(swapper, by = c('biosample_id' = 'samp')) %>% 
					mutate(bc = str_extract(NAME, '[ATGCN]+$'), Barcode = glue('{bc}_{sample_accession}'),
							Paper = 'Gautam et al. 2022') %>%
					select(Barcode, CellType, SubCellType)

SRP255012_meta <- cell_info %>% filter(study_accession == 'SRP255012') %>%
                    left_join(SRP255012, by = c('value' = 'Barcode'))

# MERGE
cell_info <- data.table::fread(config$cell_info)

meta_SRP <- bind_rows(meta_srp223254, meta_SRP158081, meta_SRP050054, meta_SRP075719, meta_MacaqueSanes, meta_SRP194595, 
						meta_mennon, meta_SRP212151, meta_mtab7316, meta_SRP257883, meta_TM, meta_SRP255195, meta_EGAD00001006350, meta_SRP218652, meta_SRP259930, 
						meta_SRP200499, outflow_meta, chick_meta, pan_human_meta, meta_SRP310237, SRP275814, SRP228556_meta, SRP255012_meta) %>%
			select(value:batch, CellType, SubCellType, TabulaMurisCellType, Paper) %>% 
	mutate(CellType = gsub('AC/HC_Precursor', 'AC/HC Precursor', CellType),
			CellType = gsub('AC/HC_Precur', 'AC/HC Precursor', CellType),
			CellType = gsub('Natural Killer', 'T/NK-Cell', CellType), 
			CellType = gsub('Artery', 'Blood Vessel', CellType),
			CellType = gsub('Vein', 'Blood Vessel', CellType),
			CellType = gsub('Erthyrocyte', 'Red Blood Cell', CellType),
			CellType = gsub('Cone Bipolar Cells', 'Bipolar Cells', CellType),
			CellType = gsub('B-cell', 'B Cell', CellType),
			CellType = gsub('T-Cell', 'T/NK-Cell', CellType),
			CellType = gsub('ium$', 'ial', CellType),
			CellType = gsub('Vascular ', '', CellType),
			CellType = gsub('s$','', CellType),
			CellType = gsub('Choriocapillari','Choriocapillaris', CellType))
cell_info_labels <- bind_rows(meta_SRP, 
                              cell_info %>% select(value:batch) %>% 
                                filter(!value %in% meta_SRP$value,
										!is.na(study_accession)) %>% 
                                mutate(Paper = NA, TissueNote = NA))


# label internal iPSC RPE data
high_ttr_huf <- read_tsv(glue('{git_dir}/data/hufnagel_iPSC_RPE_03_highTTR_barcodes.txt') )
cell_info_labels <- cell_info_labels %>% 
  mutate(CellType = case_when(#grepl('iPSC_RPE_scRNA_01', value) ~ 'iPSC',
                              #grepl('iPSC_RPE_scRNA_02', value) ~ 'RPE',
                              value %in% high_ttr_huf$value ~ 'RPE',
                              TRUE ~ CellType),
         Paper = case_when(grepl('iPSC_RPE_scRNA', value) ~ 'Swamy et al. 2021',
                           TRUE ~ Paper))


if ( sum(is.na(cell_info_labels$study_accession)) > 0 ) {
	print('Join error, check data frame!')
	save(file = 'cell_info_labels_error.Rdata')
	stop()
}
cell_info <- data.table::fread('pipeline_data/cell_info/all_cell_info.tsv') %>% select(-TissueNote) %>% filter(!is.na(study_accession))
if (((cell_info_labels$value %>% duplicated) %>% sum() == 0) & 
	(nrow(cell_info) == nrow(cell_info_labels))) {
	save(cell_info_labels, file = 'pipeline_data/cell_info/cell_info_labelled.Rdata')
	print('SUCCESS')
} else {
	save(meta_SRP218652x, meta_SRP218652y, cell_info_labels, cell_info, file = 'pipeline_data/cell_info/cell_info_labelled_fail.Rdata')
	print("Doubled or missing cells! Check data frame!")
	stop()
}


