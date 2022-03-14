# 1. hand moves over some common labels from tabula muris
hand_fixer <- function(umap){
	
	umapN <- umap #%>% 
				#mutate(CellType_predict = case_when(grepl('E12', biosample_title) & CellType_predict == 'Horizontal Cell' & organism == 'Gallus gallus' & CellType_predict_max_prob > 0.999 ~ 'Horizontal Cell',
				#								grepl('E12', biosample_title) & CellType_predict == 'Horizontal Cell' & organism == 'Gallus gallus' & CellType_predict_max_prob < 0.999 ~ 'Amacrine Cell',
				#								TRUE ~ CellType_predict)) 
	if ('TabulaMurisCellType_predict'  %in% colnames(umapN)){
	umapN <- umapN %>%
        	    mutate(CellType_predict = case_when(
									TabulaMurisCellType_predict == 'T cell' ~ 'T/NK-Cell',
                                    TabulaMurisCellType_predict == 'B cell' ~ 'B-Cell',
                                    TabulaMurisCellType_predict == 'endothelial cell' ~ 'Endothelial',
                                    TabulaMurisCellType_predict == 'epithelial cell' ~ 'Epithelial',
                                    TabulaMurisCellType_predict == 'endothelial cell' ~ 'Epithelial',
                                    TabulaMurisCellType_predict == 'keratinocyte' ~ 'Keratinocyte',
                                    TabulaMurisCellType_predict == 'blood cell' ~ 'Red Blood Cell',
                                    TabulaMurisCellType_predict == 'hepatocyte' ~ 'Hepatocyte',
                                    TabulaMurisCellType_predict == 'mesenchymal cell' ~ 'Mesenchymal',
                                    TabulaMurisCellType_predict == 'bladder cell' ~ 'Bladder',
                                    TabulaMurisCellType_predict == 'mesenchymal stem cell' ~ 'Mesenchymal (Stem)',
                                    TabulaMurisCellType_predict == 'bladder urothelial cell' ~ 'Bladder Urothelial',
                                    TabulaMurisCellType_predict == 'kidney proximal straight tubule epithelial cell' ~ 'Kidney Proximal Tubule',
                                    TabulaMurisCellType_predict == 'basal cell of epidermis' ~ 'Basal Cell',
                                    TabulaMurisCellType_predict == 'macrophage' ~ 'Macrophage',
                                    TabulaMurisCellType_predict == 'natural killer cell' ~ 'T/NK-Cell',
                                    TabulaMurisCellType_predict == 'monocyte' ~ 'Monocyte',
                                    TRUE ~ CellType_predict),
      					 CellType = case_when(TabulaMurisCellType == 'T cell' ~ 'T/NK-Cell',
                            TabulaMurisCellType == 'B cell' ~ 'B-Cell',
                            TabulaMurisCellType == 'endothelial cell' ~ 'Endothelial',
                            TabulaMurisCellType == 'epithelial cell' ~ 'Epithelial',
                            TabulaMurisCellType == 'endothelial cell' ~ 'Epithelial',
                            TabulaMurisCellType == 'keratinocyte' ~ 'Keratinocyte',
                            TabulaMurisCellType == 'blood cell' ~ 'Red Blood Cell',
                            TabulaMurisCellType == 'hepatocyte' ~ 'Hepatocyte',
                            TabulaMurisCellType == 'mesenchymal cell' ~ 'Mesenchymal',
                            TabulaMurisCellType == 'bladder cell' ~ 'Bladder',
                            TabulaMurisCellType == 'mesenchymal stem cell' ~ 'Mesenchymal (Stem)',
                            TabulaMurisCellType == 'bladder urothelial cell' ~ 'Bladder Urothelial',
                            TabulaMurisCellType == 'kidney proximal straight tubule epithelial cell' ~ 'Kidney Proximal Tubule',
                            TabulaMurisCellType == 'basal cell of epidermis' ~ 'Basal Cell',
                            TabulaMurisCellType == 'macrophage' ~ 'Macrophage',
                            TabulaMurisCellType == 'natural killer cell' ~ 'T/NK-Cell',
                            TabulaMurisCellType == 'monocyte' ~ 'Monocyte',
                            TRUE ~ CellType))
	}
	
	umapN


}
