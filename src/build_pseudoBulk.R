library(Seurat)
library(tidyverse)
library(BiocParallel)
library(SingleCellExperiment)
library(DESeq2)
library(BiocParallel)

git_dir = Sys.getenv('SCIAD_GIT_DIR')

args = commandArgs(trailingOnly=TRUE)

load(args[1]) # seurat obj
load(args[2]) # cluster
load(args[3]) # cell type prediction
group = args[4]
organism = gsub('_', ' ', args[5])
register(MulticoreParam(as.integer(args[6])))
output = args[7]

# remove non-tissue from diff testing
umap <- umap %>% filter(Source == 'Tissue')
# hand fix some labels
source(glue::glue('{git_dir}/src/tweak_celltype_labels.R'))
umap <- hand_fixer(umap)

# cut down to tissue data
integrated_obj <- integrated_obj[,umap$Barcode]
# update metadata
f_meta <- integrated_obj@meta.data  %>% as_tibble(rownames = 'Barcode') %>% select(Barcode) %>% left_join(umap) %>% left_join(meta) %>% data.frame()
row.names(f_meta) <- f_meta$Barcode
integrated_obj@meta.data <- f_meta

int_sce <-  as.SingleCellExperiment(integrated_obj)
colData(int_sce)$cluster <- as.factor(colData(int_sce)$cluster)

rm(integrated_obj)

# species level matrices
summer <- function(sce, against, species){
	print(against)
	sce_org <- sce[, umap %>% filter(organism == species) %>% pull(Barcode)]
	# remove groupings "against" that have fewer than 50 cells
	bc_to_keep <- colData(sce_org) %>% as_tibble() %>% 
		group_by(across(all_of(against))) %>% summarise(Count = n(), bc = list(Barcode)) %>% 
		filter(Count >= 50)  %>% 
		pull(bc) %>% 
		unlist()
	sce_org <- sce_org[, bc_to_keep]
	summed <- scater::aggregateAcrossCells(sce_org,
										     ids=colData(sce_org)[,against])
	summed
}

org_mat <- summer(int_sce, c(group, 'study_accession'), organism)
colData(org_mat)$cluster <- as.factor(colData(org_mat)$cluster)

deseq2_runner <- function(summed, group, formula){

	mat <- assay(summed, 'counts')
	colData(summed)$group <- colData(summed)[,group]
	colnames(mat) <-colData(summed) %>% as_tibble() %>% 
		  mutate(names = glue::glue("{study_accession}_{group}")) %>% 
			     
			pull(names) 
	# remove 0 count genes
	genes_keep <- rowSums(mat) %>% enframe() %>% filter(value > 0) %>% pull(name)
	mat <- mat[genes_keep, ]
	dds <- DESeqDataSetFromMatrix(countData = round(mat),
                              colData = colData(summed),
                              design= as.formula(formula)) # study as covariate, testing celltype diff
	# https://support.bioconductor.org/p/105087/
	dds <- DESeq(dds, parallel = TRUE, betaPrior=TRUE)
	dds
}

#deseq2_obj <- deseq2_runner(org_mat, paste0("~ 0 + ", group, " + study_accession"))
deseq2_obj <- deseq2_runner(org_mat, group, paste0("~ ", group, " + study_accession"))

save(deseq2_obj, org_mat, file = output)

#deseq_res <- results(homo_deseq2, 
#                     contrast = c("CellType_predict", "Rod", "Cone"),
#                     parallel = TRUE)
