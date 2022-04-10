# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = cbmc.adt)

# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay



library(SingleCellExperiment)
library(Seurat)
library(dplyr)
scEiaD_velocity_seurat <- CreateSeuratObject(counts = assay(velo_sce, 'spliced'), assay = 'splicedRNA')

for (i in c('unspliced','velocity','variance_velocity')){
	print(i)
	if (i %in% names(scEiaD_velocity_seurat@assays)){
		next()
	} else {
		temp <- CreateAssayObject(counts = assay(velo_sce, i))
		scEiaD_velocity_seurat[[i]] <- temp
	}
}

scEiaD_velocity_seurat@meta.data <- colData(velo_sce) %>% data.frame()



