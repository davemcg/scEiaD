args <- commandArgs(trailingOnly = TRUE)

# transfer clark ... blackshaw cell labels to all cells
library(tidyverse)
library(Seurat)


#plan(strategy = "multicore", workers = 12)
#options(future.globals.maxSize = 2400000 * 1024^2)
#load('integrated_obj.Rdata')

# load seurat integrated_obj
load(args[1])
# load umap (labelled cells)
load(args[2])
# type
transform = args[3]

nmeta <- integrated_obj@meta.data

# join with umap
nmeta <- left_join(nmeta %>% as_tibble(rownames = 'Barcode'), cell_info_labels %>% select(Barcode = value, CellType))

# label the non clark blackshaw labels as missing
nmeta$CellType[is.na(nmeta$CellType)] <- 'Missing'
integrated_obj@meta.data$CellType <- nmeta$CellType

# # var genes
if (args == 'SCT'){
  integrated_obj@assays$SCT@var.features <- grep('^MT-', integrated_obj@assays$SCT@scale.data %>% row.names(), 
                                                 value = TRUE, invert = TRUE)
} 


# # add MNN dim reduction as assay
# integrated_obj[['mnn_assay']] <- CreateAssayObject(data = t(integrated_obj@reductions$mnn@cell.embeddings))
# DefaultAssay(integrated_obj) <- 'mnn_assay'

#mnn_assay <- CreateAssayObject(data = t(integrated_obj@reductions$mnn@cell.embeddings))
#DefaultAssay(integrated_obj) <- 'integrated'
missing <- subset(integrated_obj, subset = CellType == 'Missing')
labelled <- subset(integrated_obj, subset = CellType != 'Missing')


anchors <-  FindTransferAnchors(reference = labelled, 
                                query = missing, 
                                # reference.assay = 'SCT',
                                # query.assay = 'SCT',
                                dims = 1:30)

predictions <- TransferData(anchorset = anchors, 
                            refdata = labelled@meta.data$CellType, 
                            dims = 1:30)

transferred_labels <- left_join(nmeta, 
                                predictions %>% as_tibble(rownames = 'Barcode')) %>% 
  mutate(ID = case_when(CellType == 'Missing' ~ `predicted.id`, 
                        TRUE ~ CellType))

#integrated_obj@meta.data$CellType_transfer <- meta %>% pull(ID)

# save(integrated_obj, file = 'integrated_obj__transfer.Rdata', compress = FALSE)
#save(integrated_obj, file = args[4], compress = FALSE)
save(predictions, file = args[4])