args <- commandArgs(trailingOnly = TRUE)

# transfer clark ... blackshaw cell labels to all cells
library(tidyverse)
library(Seurat)

print(args)

#plan(strategy = "multicore", workers = 12)
#options(future.globals.maxSize = 2400000 * 1024^2)
#load('integrated_obj.Rdata')

# load well seurat integrated_obj
load(args[1])
well_obj <- integrated_obj
# load droplet obj
load(args[2])
droplet_obj <- integrated_obj

# load labelled cells
load(args[3])
# type
transform = args[4]

# label cells
labeller <- function(seurat_obj){
	nmeta <- seurat_obj@meta.data
	# join with cell info
	nmeta <- left_join(nmeta %>% as_tibble(rownames = 'Barcode'), cell_info_labels %>% select(Barcode = value, CellType))
	# label the unlabelled cells as missing
	nmeta$CellType[is.na(nmeta$CellType)] <- 'Missing'
	seurat_obj@meta.data$CellType <- nmeta$CellType
	seurat_obj
}

well_obj <- labeller(well_obj)
droplet_obj <- labeller(droplet_obj)


# check var genes and swap SCT to RNA
seurat_fixer <- function(seurat_obj){
	# # var genes
	if (seurat_obj@assays[[DefaultAssay(seurat_obj)]]@var.features %>% length() == 0){
 	 seurat_obj@assays[[DefaultAssay(seurat_obj)]]@var.features <- 
   	 grep('^MT-', 
   	      seurat_obj@assays[[DefaultAssay(seurat_obj)]]@scale.data %>% row.names(), 
   	      value = TRUE, invert = TRUE)
	} 

	# if assay is SCT there's some kind of bug in RenameCells
	# so for now I'm just going to swap over to RNA
	if (DefaultAssay(seurat_obj) == 'SCT' || DefaultAssay(seurat_obj) == 'integrated') {
		seurat_obj@assays$RNA@var.features = seurat_obj@assays[[DefaultAssay(seurat_obj)]]@var.features 
		DefaultAssay(seurat_obj) <- 'RNA'
		seurat_obj@assays$SCT <- NULL
		for (i in seq(1,length(seurat_obj@reductions))){
 			 seurat_obj@reductions[[i]] <- NULL
		}
		seurat_obj@reductions$SCT <- NULL
		seurat_obj@reductions[[1]] <- NULL
	}
	seurat_obj
}

well_obj <- seurat_fixer(well_obj)
droplet_obj <- seurat_fixer(droplet_obj)

droplet_obj <- subset(droplet_obj, subset = CellType != 'Missing')

# further subset droplet obj to ~50k cells
# no more than 4000 cells per CellType
# was getting memory errors
set.seed(421)
to_keep <- droplet_obj@meta.data  %>% as_tibble(rownames = 'cells') %>% 
	filter(!CellType %in% c('Doublet', 'Doublets')) %>%
	group_by(CellType) %>% 
	sample_n(5000, replace = TRUE) %>% 
	unique()  %>% 
	pull(cells)
droplet_obj <- subset(droplet_obj, cells = to_keep)


anchors <-  FindTransferAnchors(reference = droplet_obj, 
                                query = well_obj) 
                                # reference.assay = 'SCT',
                                # query.assay = 'SCT',
                                #dims = 1:30)

predictions <- TransferData(anchorset = anchors, 
                            refdata = droplet_obj@meta.data$CellType)
                            #dims = 1:30)

save(predictions, file = args[5])
