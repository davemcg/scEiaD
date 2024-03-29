# extract UMAP coords labelled with meta data
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(tidyverse)
# read in ENS <-> HGNC gene mapping info
gene_map <- read_tsv('references/ENSG2gene_name.tsv.gz') %>% mutate(hs_gene_name = toupper(hs_gene_name))

# load cluster data
load(args[3])
cluster <- meta %>% pull(2)
subcluster <- meta %>% pull(3)
# load pre int seurat obj, calc cell cycle
load(args[1])
# convert Seurat cell cycle HGNC to ENSGENE
s.genes <- cc.genes$s.genes %>% enframe(value = 'hs_gene_name') %>% left_join(gene_map) %>% pull(hs_gene_id) 
g2m.genes <- cc.genes$g2m.genes %>% enframe(value = 'hs_gene_name') %>% left_join(gene_map) %>% pull(hs_gene_id)
if (DefaultAssay(integrated_obj) == 'integrated'){
	DefaultAssay(integrated_obj) <- 'RNA'
}
integrated_obj <- CellCycleScoring(integrated_obj, s.features = s.genes, g2m.features = g2m.genes)
meta <- integrated_obj@meta.data
# load integrated seurat obj
load(args[2])
integrated_obj@meta.data$`S.Score` <- meta$`S.Score`
integrated_obj@meta.data$`G2M.Score` <- meta$`G2M.Score`
integrated_obj@meta.data$`Phase` <- meta$`Phase`
integrated_obj@meta.data$cluster <- cluster
integrated_obj@meta.data$subcluster <- subcluster
# load labelled cell data
load(args[4])
# load predicted cell data (+ labelled)
#load(args[4])
# method
method <- args[6]

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

if (args[6] == 'TSNE'){
  reduction.key <- gsub('UMAP','TSNE', reduction.key)
  reduction.name <- gsub('UMAP','TSNE', reduction.name)
}

# left_join known cell labels
orig_meta <- integrated_obj@meta.data %>% as_tibble(rownames = 'Barcode')
umap <- Embeddings(integrated_obj[[reduction.name]]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  left_join(., cell_info_labels %>% select(-contains(c('study_accession', 'Age', 'batch'))) %>% rename(Barcode = value),
            by = 'Barcode') 
#  left_join(., predictions %>% 
#              as_tibble(rownames = 'Barcode') %>% 
#              select(Barcode, CellType_predict = `predicted.id`)) %>% 
#  mutate(CellType_predict = case_when(is.na(CellType_predict) ~ CellType,
#                                      TRUE ~ CellType_predict))

umap$Method <- method
if (args[7] == 'UMAP'){
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
#DefaultAssay(integrated_obj) <- 'RNA'
#core_expression <- FetchData(integrated_obj, core_markers, slot = 'counts') %>% as_tibble(rownames = 'Barcode')

#umap <-   left_join(umap, core_expression, by = 'Barcode')

save(umap, file = args[5])
