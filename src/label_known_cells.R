library(tidyverse)
load('~/git/massive_integrated_eye_scRNA/data/sra_metadata_extended.Rdata')
meta <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv') %>% select(-TissueNote)
# load labelled data from clark et all
# NO WRONG NOW https://www.dropbox.com/s/y5lho9ifzoktjcs/10x_mouse_retina_development_phenotype.csv?dl=1
clark_labels <- data.table::fread('GSE118614_barcodes.tsv.gz')
#clark_labels <- read_csv('10x_mouse_retina_development_phenotype.csv')

# extract clark blackshaw fields we1 want
clark_labels <- clark_labels %>% 
  mutate(UMI = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>% 
  select(CellType = umap2_CellType, umap_coord1, 
         umap_coord2, umap_coord3, umap_cluster, V1, barcode, sample, age, 
         UMI)

## now get macosko labels
# macosko et al
macosko_labels <- read_tsv('http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt', col_names = c('Cell','Cluster'))
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
cell_info <- read_tsv('cell_info.tsv') %>% select(-TissueNote)
cell_info <- cell_info %>% filter(study_accession == 'E-MTAB-7316') %>%
				mutate(UMI = gsub('_\\w+', '', value)) %>%
				mutate(sample = case_when(sample_accession == 'ERS2852885' ~ '5',
											sample_accession == 'ERS2852886' ~ '3',
											sample_accession == 'ERS2852887' ~ '4',
											sample_accession == 'ERS2852888' ~ '1',
											sample_accession == 'ERS2852889' ~ '2'))
												
retina_wong <- read_csv('~/git/massive_integrated_eye_scRNA/data/retina_wong_cellbc_cellid.csv') %>%
				mutate(CellType = gsub(' C\\d+', '', cell.id.orig)) %>% dplyr::select(`cell.bc`, CellType) %>%
				mutate(CellType = case_when(grepl('RGC', CellType) ~ 'Retinal Ganglion Cells',
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
rgc_crush <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sanes__RGC_Atlas_coordinates.txt')
rgc_crush <- rgc_crush[-1,]
rgc_crush <- rgc_crush %>% 
				separate(NAME, c('sample', 'UMI'), sep = '_') %>%
				mutate(UMI = gsub('-\\d+','',UMI))
rgc_crush$CellType <- 'Retinal Ganglion Cells'
## load cell info
cell_info <- read_tsv('cell_info.tsv') %>% select(-TissueNote)
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
karthik <- read_tsv('~/git/massive_integrated_eye_scRNA/data/shekhar_sanes_ClustAssignFile.txt')
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
cell_info <- read_tsv('cell_info.tsv') %>% select(-TissueNote)
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
cell_info <- read_tsv('cell_info.tsv')  
lu_clark <- data.table::fread('~/git/massive_integrated_eye_scRNA/data/GSE138002_Final_barcodes.csv.gz') %>%
							mutate(UMI = gsub('^.*\\.','', barcode) %>% gsub('-.*','',.)) %>%
				dplyr::rename(CellType = umap2_CellType) %>%
				mutate(CellType = gsub('BC/Photo_Precurs', 'Photoreceptor Precursors', CellType))
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
  left_join(lu_clark %>% select(UMI, sample, CellType),
				by = c('UMI', 'sample')) %>% 
  select(value:batch,CellType) %>% mutate(Paper = 'Lu et al. 2020')
## sanes macaque 
cell_info <- read_tsv('cell_info.tsv') %>% select(-TissueNote)
sanes_files <- list.files('~/git/massive_integrated_eye_scRNA/data/', "Macaque*", full.names = TRUE)
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
scheetz_files <- list.files('~/git/massive_integrated_eye_scRNA/data/', pattern = 'GSM.*donor.*gz', full.names = TRUE)
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
  mutate(Paper = 'Voigt et al. 2019')

# mennon et al. 
# SRP222001, SRP222958
mennon_seqwell <- read_tsv('~/git/massive_integrated_eye_scRNA/data/GSE137846_Seq-Well_sample_annotations.txt.gz', 
                           col_names = FALSE, skip = 1) %>% 
  dplyr::select(barcode = X1, Labels = X49)
mennon_drop <- read_tsv('~/git/massive_integrated_eye_scRNA/data/GSE137537_sample_annotations.tsv.gz') %>% 
  dplyr::select(barcode = Barcode, tissue, Labels)
mennon <- bind_rows(mennon_seqwell, mennon_drop) %>% 
  mutate(CellType = case_when(Labels == 'ACs' ~ 'Amacrine Cells',
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
load('~/git/massive_integrated_eye_scRNA/data/SRP257883_GSE_meta.Rdata')
SRP257883_GSE_meta <- SRP257883_GSE_meta %>% 
  mutate(CellType = case_when(celltype == 'smc' ~ 'Smooth Muscle Cell',
                              celltype == 'RBC' ~ 'Red Blood Cells',
                              celltype == 'fibroblast' ~ 'Fibroblasts',
                              celltype == 'macrophage' ~ 'Macrophages',
                              TRUE ~ stringr::str_to_title(celltype))) %>% 
  filter(!CellType %in% c('Rod_rpe')) %>% 
  select(-celltype)
cell_info <- read_tsv('cell_info.tsv') %>%  
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
# SRP218652 mullins/scheetz RPE
cell_info <- read_tsv('cell_info.tsv') %>%
    filter(study_accession == 'SRP218652') %>%
    mutate(UMI = gsub('_\\w+', '', value))
load('~/git/massive_integrated_eye_scRNA/data/SRP218652__meta.Rdata')
SRP218652 <- SRP218652 %>%  select(Barcode, sample_accession, CellType) %>% 
	mutate(UMI =  gsub('_\\w+', '', Barcode),
			UMI = gsub('-\\d+', '', UMI))
meta_SRP218652 <- cell_info %>% left_join(., SRP218652, by = c('UMI', 'sample_accession'))

# tabula muris
cell_info <- read_tsv('cell_info.tsv') %>% select(-TissueNote) %>% filter(study_accession == 'SRP131661') 
tm <- data.table::fread('~/git/massive_integrated_eye_scRNA/data/tabula_muris_meta.tsv.gz') %>%
			dplyr::rename(TabulaMurisCellType = CellType) %>% 
  			mutate(CellType = case_when(TabulaMurisCellType == 'B cell' ~ 'B-Cell',
                              TabulaMurisCellType == 'blood cell' ~ 'Red Blood Cells',
                              TabulaMurisCellType == 'T cell' ~ 'T-Cell',
                              TabulaMurisCellType == 'fibroblast' ~ 'Fibroblasts',
                              grepl('endoth', TabulaMurisCellType, ignore.case = TRUE) ~ 'Endothelial')) %>% 
			select(value, Paper, CellType, TabulaMurisCellType)
meta_TM <- cell_info %>% left_join(., tm, by = 'value' )

## Yan ... Sanes retina rod/cone 
cell_info <- read_tsv('cell_info.tsv') %>%
    filter(study_accession == 'SRP255195')
sample_meta <- read_tsv('~/git/massive_integrated_eye_scRNA/data/sample_run_layout_organism_tech.tsv') %>% filter(study_accession == 'SRP255195')
yan_pr <- read_csv('~/git/massive_integrated_eye_scRNA/data/yan_sanes_PR____Human_retina_combined_all_meta.csv')
yan_pr <- yan_pr %>% 
  mutate(CellType = case_when(Cluster == 'Astrocytes' ~ 'Astrocytes',
                              grepl('^DB|FMB|IMB|RB1|OFFx|BB', Cluster) ~ 'Bipolar Cells',
                              Cluster == 'Endothelium' ~ 'Endothelium',
                              grepl('Gaba|Gly', Cluster) ~ 'Amacrine Cells',
                              grepl('H1|H2', Cluster) ~ 'Horizontal Cells',
                              grepl('MG_|Muller', Cluster) ~ 'Muller Glia',
                              grepl('MG_OFF|MG_ON|RGC|PG_|IMB', Cluster) ~ 'Retinal Ganglion Cells',
                              grepl('Cone', Cluster) ~ 'Cone',
                              grepl('Rod', Cluster) ~ 'Rods',
                              Cluster == 'MicroGlia' ~ 'Microglia'),
         SubCellType = Cluster,
         sample = str_extract(NAME, 'H\\d+\\w+S\\d')) %>% 
  filter(!is.na(CellType)) %>% left_join(sample_meta %>%
                                           mutate(sample = str_extract(TissueNote, 'H\\d+\\w+S\\d')), 
                                         by = 'sample') %>% 
  filter(!is.na(sample_accession)) %>% 
  separate(NAME, into = c('samp','UMI','lane'), sep = "[_|-]") %>% 
  mutate(value = paste0(sample_accession, '_', UMI)) %>%
  mutate(Paper = 'Yan et al. 2020')
meta_SRP255195 <- cell_info %>% left_join(., yan_pr %>% select(value, Paper, CellType, SubCellType), by = 'value')


meta_SRP <- bind_rows(meta_SRP218652, meta_srp223254, meta_SRP158081, meta_SRP050054, meta_SRP075719, meta_MacaqueSanes, meta_SRP194595, 
						meta_mennon, meta_SRP212151, meta_mtab7316, meta_SRP257883, meta_TM, meta_SRP255195) %>%
	mutate(CellType = gsub('Melanotype', 'Melanocytes', CellType),
			CellType = gsub('B-cell', 'B-Cell', CellType),
			CellType = gsub('Macrophages', 'Macrophage', CellType),
			CellType = gsub('Pericyte$', 'Pericytes', CellType),
			CellType = gsub('T/NK-cell', 'T-Cell', CellType))
cell_info_labels <- bind_rows(meta_SRP, 
                              cell_info %>% select(value:batch) %>% 
                                filter(!value %in% meta_SRP$value) %>% 
                                mutate(Paper = NA, TissueNote = NA))
# this is crude, but SRP16660 has selected Muller Glia via FACS of R26R mouse line
# SRP186407 uses Cx3cr1+ FACS to select Microglia
cell_info_labels <- cell_info_labels %>% 
  mutate(CellType = case_when(study_accession == 'SRP166660' ~ "Muller Glia",
                              study_accession == 'SRP186407' ~ "Microglia",
                              TRUE ~ CellType),
         Paper = case_when(study_accession == 'SRP166660' ~ "Rueda et al. 2019",
                           study_accession == 'SRP186407' ~ "O'Koren et al. 2019",
                           TRUE ~ Paper))


# label internal iPSC RPE data
cell_info_labels <- cell_info_labels %>% 
  mutate(CellType = case_when(grepl('iPSC_RPE_scRNA_01', value) ~ 'iPSC',
                              grepl('iPSC_RPE_scRNA_02', value) ~ 'RPE',
                              grepl('iPSC_RPE_scRNA_03', value) ~ 'RPE (Transwell)',
                              TRUE ~ CellType),
         Paper = case_when(grepl('iPSC_RPE_scRNA', value) ~ 'Hufnagel 2020',
                           TRUE ~ Paper))

if ((cell_info_labels$value %>% duplicated) %>% sum() == 0) {
	save(cell_info_labels, file = 'cell_info_labelled.Rdata')
} else {
	print("Doubled cell labels! Check data frame!")
}


# ## Figure out what the hell is going on as the above file lists p1, r1 through r6 and GEO lists the 7 samples as 1 - 7....
# ## I'm guessing either goes p1, r1, r2...r6 or r1...r6, p1. 
# ## will use barcode overlap to figure this out
# macosko_barcodes <- colnames(fread('GSE63472_P14Retina_merged_digital_expression.txt.gz', nrows = 1))
# macosko_barcodes <- macosko_barcodes[2:length(macosko_barcodes)]
# ## load "mouse retina 1" and figure out what barcode set it matches best....
# ## SRS866912
# mr1 <- row.names(seurat_late__standard@meta.data) %>% enframe() %>% filter(grepl('SRS866912', value)) %>% mutate(umi = gsub('_.*||\\.\\d*', '', value))
# ## this returns 6450
# ## so retina 1 is r1
# grep('r1_', macosko_barcodes, value = T) %>% gsub('.*_', '', .) %>% enframe() %>% filter(value %in% mr1$umi) %>% dim()
# ## let's see if retina 7 is p1
# ## it is!
# mr7 <- row.names(seurat_late__standard@meta.data) %>% enframe() %>% filter(grepl('SRS866906', value)) %>% mutate(umi = gsub('_.*||\\.\\d*', '', value))
# grep('r7_', macosko_barcodes, value = T) %>% gsub('.*_', '', .) %>% enframe() %>% filter(value %in% mr7$umi) %>% dim()
# ## ok, so appears to p1 == r7
# ## let's just check one in the middle to be sure, r4
# ## yep
# mr4 <- row.names(seurat_late__standard@meta.data) %>% enframe() %>% filter(grepl('SRS866909', value)) %>% mutate(umi = gsub('_.*||\\.\\d*', '', value))
# grep('r4_', macosko_barcodes, value = T) %>% gsub('.*_', '', .) %>% enframe() %>% filter(value %in% mr4$umi) %>% dim()
