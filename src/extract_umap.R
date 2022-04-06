args= c('/data/swamyvs/scEiaD/rson_tmp/_7dbjtmv.json', 'UMAP')
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(tidyverse)
# read in ENS <-> HGNC gene mapping info
gene_map <- read_tsv('references/ENSG2gene_name.tsv.gz') %>% mutate(hs_gene_name = toupper(hs_gene_name))
library(jsonlite)
rule <- read_json(args[1])
# read in ENS <-> HGNC gene mapping info
partition <- rule$wildcards$partition

if(partition %in% c('human',  'mouse','chick','macaque')){
  gene_map <- rtracklayer::readGFF(rule$input$gene_id_mapper) %>% 
    as_tibble %>% 
    filter(type == 'transcript') %>% 
    select(gene_id, gene_name) %>% distinct %>% 
    mutate(gene_name = toupper(gene_name),
           gene_id = str_split(gene_id ,'\\.') %>% sapply(function(x) x[1]))
} else {
  gene_map <- read_tsv(rule$input$gene_id_mapper) %>% 
    mutate(hs_gene_name = toupper(hs_gene_name))  %>% 
    rename(gene_id = hs_gene_id, gene_name = hs_gene_name)
}


# load cluster data
load(rule$input$cluster_rdata)
cluster <- meta %>% pull(2)
subcluster <- meta %>% pull(3)
# load pre int seurat obj, calc cell cycle
# load int seurat obj, calc cell cycle
load(rule$input$intg_seu_obj)
# convert Seurat cell cycle HGNC to ENSGENE
s.genes <- cc.genes$s.genes %>% enframe(value = 'gene_name') %>% left_join(gene_map) %>% pull(gene_id) 
g2m.genes <- cc.genes$g2m.genes %>% enframe(value = 'gene_name') %>% left_join(gene_map) %>% pull(gene_id)
if (DefaultAssay(integrated_obj) == 'integrated'){
  DefaultAssay(integrated_obj) <- 'RNA'
}
integrated_obj <- CellCycleScoring(integrated_obj, s.features = s.genes, g2m.features = g2m.genes)
meta <- integrated_obj@meta.data
# load umap seurat obj.
load(rule$input$umap_seu_obj)
# remove dup barcode cols
bc_cols <- grep
integrated_obj@meta.data$`S.Score` <- meta$`S.Score`
integrated_obj@meta.data$`G2M.Score` <- meta$`G2M.Score`
integrated_obj@meta.data$`Phase` <- meta$`Phase`
integrated_obj@meta.data$cluster <- cluster
integrated_obj@meta.data$subcluster <- subcluster
# load labelled cell data
load(rule$input$cell_info_labeled)
# load predicted cell data (+ labelled)
#load(args[4])
# method
method <- rule$wildcards$method

if (method == 'CCA'){
  reduction <- 'pca'
  reduction.key <- 'ccaUMAP_'
  #DefaultAssay(integrated_obj) <- 'SCT'
} else if (method == 'fastMNN'){
  reduction <- 'mnn'
  reduction.key <- 'mnnUMAP_'
} else if (method == 'none'){
  reduction <- 'pca'
  reduction.key <- 'noneUMAP_'
} else if (method == 'combat'){
  if ((integrated_obj@reductions$pca@cell.embeddings %>% ncol()) < 100){
    integrated_obj <- RunPCA(integrated_obj, npcs = 100)
  }
  reduction <- 'pca'
  reduction.key <- 'combatUMAP_'
} else if (method == 'liger'){
  reduction <- 'iNMF'
  reduction.key <- 'iNMFUMAP_'
} else if (grepl('scVI|scANV', method)) {
   reduction <- 'scVI'
   reduction.key <- 'scviUMAP_'
} else {
  print("GUESSING!")
  reduction <- method
  reduction.key <- paste0(reduction, 'UMAP_')
}
reduction.name <- gsub('_','', reduction.key)

if (rule$wildcards$method == 'TSNE'){
  reduction.key <- gsub('UMAP','TSNE', reduction.key)
  reduction.name <- gsub('UMAP','TSNE', reduction.name)
}

# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(integrated_obj[[reduction.name]]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., cell_info_labels %>% select(-contains(c('study_accession', 'Age', 'batch'))) %>% rename(Barcode = value),
            by = 'Barcode') %>%
  mutate(CellType = gsub('AC/HC_Precursorsor', 'AC/HC_Precursor', CellType))
#  left_join(., predictions %>% 
#              as_tibble(rownames = 'Barcode') %>% 
#              select(Barcode, CellType_predict = `predicted.id`)) %>% 
#  mutate(CellType_predict = case_when(is.na(CellType_predict) ~ CellType,
#                                      TRUE ~ CellType_predict))

umap$Method <- method
if (args[2] == 'UMAP'){
  colnames(umap)[2:3] <- c('UMAP_1', 'UMAP_2')
} else {
  colnames(umap)[2:3] <- c('TSNE_1', 'TSNE_2')
}

# https://www.sciencedirect.com/science/article/pii/S1534580720303075?via%3Dihub
# In mice, the transition from early and late-stage RPCs occurs rapidly between E16 and E18 (Clark et al., 2019), but in humans this process occurs between 11 and 15 GW (Figure S2B) and likely reflects the differences in the timing of human retina development between central and peripheral regions (Diaz-Araya and Provis, 1992, van Driel et al., 1990).
umap <- umap %>% mutate(integration_group = case_when(organism == 'Homo sapiens' & Age <= -175 ~ 'Early Dev.',
                                                      organism == 'Homo sapiens' & Age <= 0 ~ 'Late Dev.',
                                                      organism == 'Homo sapiens' & Age <= 360 ~ 'Maturing',
                                                      organism == 'Homo sapiens' ~ 'Mature', 
                                                      organism == 'Mus musculus' & Age < -2 ~ 'Early Dev.',
                                                      organism == 'Mus musculus' & Age <= 0 ~ 'Late Dev.',
                                                      organism == 'Mus musculus' & Age < 14 ~ 'Maturing',
                                                      organism == 'Mus musculus' ~ 'Mature',
                                                      organism == 'Macaca fascicularis' ~ 'Mature'))

# add core markers
# Rho for rods
# Opn1sw for cones
# Sfrp2 for early progenitors
# Hes6 for neurogenic
# Tfap2a for amacrine
# Ccnd1 for late progenitors
# Aqp4 for muller glia
# Vsx1, Grm6, Cabp5 for bipolar cells
# Elavl4 for RGC
# CX3CR1 for microglia

# https://www.biorxiv.org/content/biorxiv/early/2019/07/16/703348.full.pdf
# rods: GNGT1, NRL, PDE6G, RHO; 
# cones: ARR3, CNGA3, OPN1LW, OPN1MW, OPN1SW, PDE6H; 
# horizontal cells: LHX1, LNP1, ONECUT1, VAT1L; 
# bipolar cells: GRIK1, IRX6, LRTM1, PCP2, PRKCA, TRPM1, VSX1, VSX2; 
# amacrine cells: GAD1, SLC6A9, TFAP2A, TFAP2B; 
# ganglion cells: POU4F2, NEFL, NEFM, RBPMS, SLC17A6, SNCG, THY1; 
# pigmented cells: BEST1, MITF, MLANA, TJP1, RPE65; 
# glial cells: CRABP1, GFAP, GLUL; 
# endothelial cells, mural cells and fibroblasts: ACTA2, COL1A2, EGFL7, PDGFRB, PROCR, VWF; 
# immune cells: AIF1, CD2, CD48, CX3CR1, HBB, IL32, JCHAIN, LST1.
#core_markers <- c(c('Rho','Opn1sw', 'Sfrp2', 'Hes6', 'Tfap2a', 'Isl1', 'Ccnd1','Aqp4', 'Vsx1', 'Elavl4', 'Best1', 'Grm6', 'Cabp5', "Cx3cr1"),
#                  c('GNGT1','NRL','PDE6G','RHO','ARR3','CNGA3','OPN1LW','OPN1MW','OPN1SW','PDE6H','LHX1','LNP1','ONECUT1','VAT1L','GRIK1','IRX6','LRTM1','PCP2','PRKCA','TRPM1','VSX1','VSX2','GAD1','SLC6A9','TFAP2A','TFAP2B','POU4F2','NEFL','NEFM','RBPMS','SLC17A6','SNCG','THY1','BEST1','MITF','MLANA','TJP1','RPE65','CRABP1','GFAP','GLUL','ACTA2','COL1A2','EGFL7','PDGFRB','PROCR','VWF','AIF1','CD2','CD48','CX3CR1','HBB','IL32','JCHAIN','LST1', 'PRKCA', 'SCGN', 'NTNG1', 'SNCG')) %>% 
#  toupper() %>% unique()
#core_markers <- c('ENSG00000196517','ENSG00000162374','ENSG00000116745','ENSG00000162631','ENSG00000116824','ENSG00000117091','ENSG00000135821','ENSG00000173267','ENSG00000107796','ENSG00000244734','ENSG00000091664','ENSG00000167995','ENSG00000110092','ENSG00000154096','ENSG00000110799','ENSG00000139053','ENSG00000129535','ENSG00000119614','ENSG00000104067','ENSG00000134160','ENSG00000169856','ENSG00000166426','ENSG00000008517','ENSG00000159387','ENSG00000171724','ENSG00000273706','ENSG00000131095','ENSG00000154229','ENSG00000185527','ENSG00000171885','ENSG00000174788','ENSG00000105507','ENSG00000144191','ENSG00000128683','ENSG00000144485','ENSG00000100987','ENSG00000101000','ENSG00000171189','ENSG00000168329','ENSG00000144771','ENSG00000187098','ENSG00000206535','ENSG00000163914','ENSG00000132465','ENSG00000151615','ENSG00000145423','ENSG00000016082','ENSG00000113721','ENSG00000113262','ENSG00000137203','ENSG00000079689','ENSG00000204482','ENSG00000204472','ENSG00000008196','ENSG00000127928','ENSG00000164692','ENSG00000128617','ENSG00000104722','ENSG00000277586','ENSG00000157110','ENSG00000120215','ENSG00000172889','ENSG00000277401','ENSG00000274965','ENSG00000274577','ENSG00000235915','ENSG00000226182','ENSG00000237727','ENSG00000223465','ENSG00000234836','ENSG00000231048','ENSG00000234514','ENSG00000235985','ENSG00000206433','ENSG00000206428','ENSG00000230791','ENSG00000235588','ENSG00000120500','ENSG00000102076','ENSG00000268221','ENSG00000285493')
#DefaultAssay(integrated_obj) <- 'RNA'
#core_expression <- FetchData(integrated_obj, core_markers, slot = 'counts') %>% data.frame() %>% as_tibble(rownames = 'Barcode')

#umap <-   left_join(umap, core_expression, by = 'Barcode')

save(umap, file = rule$output$umap_data )

