library(pool)
library(RSQLite)
library(Seurat)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

git_dir = Sys.getenv('SCIAD_GIT_DIR')

load(args[1])
load(args[2])
load(args[3])
load(args[4])
load(args[5])

# add study meta
study_meta <- read_tsv(paste0(git_dir, '/data/GEO_Study_Level_Metadata.tsv'))
umap <- umap %>% left_join(phate_2D$embedding %>% as_tibble(rownames = 'Barcode'))
umap <- umap %>% left_join(., study_meta, by = c('study_accession'))
well_metadata <- seurat_obj@meta.data %>% as_tibble() %>% select(-contains('RNA_snn'), -barcode, -cluster) %>% rename(cluster = seurat_clusters, CellType_predict = CellType)
well_metadata <- well_metadata %>% left_join(., study_meta, by = c('study_accession'))

# make metadata
metadata <- umap %>% select(Barcode, UMAP_1, UMAP_2, PHATE1, PHATE2, nCount_RNA, nFeature_RNA, Phase, percent_mt = `percent.mt`, S_Score = `S.Score`, G2M_Score = `G2M.Score`, batch, sample_accession, study_accession, Age, library_layout, organism, Platform, UMI, Covariate, CellType, CellType_predict, TabulaMurisCellType, TabulaMurisCellType_predict, SubCellType, GSE, Summary, Citation, PMID, Design) %>% 
		mutate(SubCellType = gsub('p_','', SubCellType)) %>% mutate(SubCellType = gsub('f_','', SubCellType)) %>%
			left_join(., meta, by = 'Barcode')

metadata$SubCellType[grepl('^RB|Rods|Peri|^MG$|^Mic$', metadata$SubCellType)] <- NA
metadata$CellType[grepl('Doub|\\/Margin\\/Periocular', metadata$CellType)] <- NA


colnames(metadata)[ncol(metadata)] <- 'subcluster'
colnames(metadata)[(ncol(metadata) - 1)] <- 'cluster'

metadata <- bind_rows(metadata , well_metadata %>% 
		select(Barcode, UMAP_1, UMAP_2, nCount_RNA, nFeature_RNA, Phase, percent_mt = `percent.mt`, S_Score = `S.Score`, G2M_Score = `G2M.Score`, batch, sample_accession, study_accession, Age, library_layout, organism, Platform, UMI, Covariate, CellType_predict, SubCellType, GSE, Summary, Citation, PMID, Design, cluster) %>% 
		mutate(cluster = as.numeric(as.character(cluster) )))  %>%  
  mutate(cluster = as.character(cluster),
			PMID = as.character(PMID)) %>%
  mutate(Stage = case_when(organism == 'Homo sapiens' & Age <= -175 ~ 'Early Dev.',
                                                organism == 'Homo sapiens' & Age <= 0 ~ 'Late Dev.',
                                                organism == 'Homo sapiens' & Age <= 360 ~ 'Maturing',
                                                organism == 'Homo sapiens' ~ 'Mature',
                                                organism == 'Mus musculus' & Age < -2 ~ 'Early Dev.',
                                                organism == 'Mus musculus' & Age <= 0 ~ 'Late Dev.',
                                                organism == 'Mus musculus' & Age < 14 ~ 'Maturing',
                                                organism == 'Mus musculus' ~ 'Mature',
                                                organism == 'Macaca fascicularis' ~ 'Mature'))

meta_filter <- metadata %>% 
  filter(!is.na(study_accession), 
         !CellType_predict %in% c('Doublet', 'Doublets'))

cpm_drop <- RelativeCounts(integrated_obj@assays$RNA@counts, scale.factor= 1e6)
cpm_well <- RelativeCounts(seurat_obj@assays$RNA@counts, scale.factor= 1e6)

if (row.names(cpm_drop) == row.names(cpm_well)) {
    print('merging droplet and well cpm')
	cpm = cbind(cpm_drop, cpm_well)
} else {
    stop('Row names do not line up!')
}

chunk_num = 20

print(dim(cpm))		
genes <- row.names(cpm) %>% enframe(value = 'Gene') %>% dplyr::select(-name) %>% arrange()

long_data <- list()

if (args[7] == 'TRUE'){
	seurat_meta <- integrated_obj@meta.data
	seurat_meta <- seurat_meta %>% as_tibble(rownames = 'Barcode') %>% left_join(., umap %>% select(Barcode, organism, CellType_predict), by = 'Barcode')
	for (i in seq(1:chunk_num)){
		chunk <- round(nrow(cpm)/chunk_num)
   		end = chunk * i
    	start = (chunk * i) - chunk + 1
    	if (i == chunk_num){
        	end = nrow(cpm)
    	}

    	print(paste(start, end))
		design = model.matrix(~ as.factor(seurat_meta$batch) + as.factor(seurat_meta$CellType_predict))
								#as.numeric(seurat_meta$`percent.mt`) + 
								#as.numeric(seurat_meta$`nCount_RNA`) + 
								#as.factor(seurat_meta$organism) +
								#as.factor(seurat_meta$CellType_predict))
    	cor_mat <- limma::removeBatchEffect(VISION:::matLog2(cpm[start:end,]), 
							#batch = as.factor(seurat_meta$batch),
							design = design)
		long_data[[i]] <- cor_mat %>%
    		as_tibble(rownames = 'Gene') %>%
    		pivot_longer(cols = 2:(ncol(cor_mat) + 1), names_to = 'Barcode', values_to = 'cpm') %>%
    		filter(cpm > 2.5)
		rm(cor_mat)
		}
	rm(seurat_meta)
} else {
	# cpm <- VISION:::matLog2(cpm)
	for (i in seq(1:chunk_num)){
		chunk <- round(ncol(cpm)/chunk_num)
		end = chunk * i
 		start = (chunk * i) - chunk + 1
 		if (i == chunk_num){
			end = ncol(cpm)
		}
		
		print(paste(start, end))
		long_data[[i]] <- cpm[,start:end] %>%
			as.matrix() %>%
			as_tibble(rownames = 'Gene') %>%
			pivot_longer(cols = 2:(end - start + 2), names_to = 'Barcode', values_to = 'cpm') %>%
			filter(cpm > 0) %>% 
			mutate(cpm = log2(cpm + 1))
	}
}
rm(integrated_obj)

long <- bind_rows(long_data)
pool <- dbPool(RSQLite::SQLite(), dbname = args[6])
dbWriteTable(pool, "cpm", long, overwrite = TRUE)
db_create_index(pool, table = 'cpm', columns = c('Gene'))
db_create_index(pool, table = 'cpm', columns = c('Barcode'))




dbWriteTable(pool, "metadata", metadata, overwrite = TRUE)
db_create_index(pool, table = 'metadata', columns = c('Barcode'))

dbWriteTable(pool, "metadata_filter", meta_filter, overwrite = TRUE)
db_create_index(pool, table = 'metadata_filter', columns = c('Barcode'))

dbWriteTable(pool, "genes", genes, overwrite = TRUE)


long <- left_join(long, meta_filter %>% select(Barcode, batch, study_accession, library_layout, organism,Stage, Platform, Covariate, CellType, SubCellType, CellType_predict, TabulaMurisCellType, TabulaMurisCellType_predict, Phase, GSE, Citation, PMID, cluster), by = 'Barcode') 

grouped_stats <- long %>% group_by(Gene, batch, study_accession, library_layout, organism,Stage, Platform, Covariate, CellType, SubCellType, CellType_predict, TabulaMurisCellType, TabulaMurisCellType_predict, Phase, GSE, Citation, PMID, cluster) %>%
						summarise(cell_ct = n(), cell_exp_ct = sum(cpm > 0), cpm = mean(cpm))
dbWriteTable(pool, 'grouped_stats', grouped_stats, overwrite = TRUE)
db_create_index(pool, table = 'grouped_stats', columns = c('Gene'))

#meta_only_grouped_stats <- long %>% group_by(batch, study_accession, library_layout, organism,Stage, Platform, Covariate, CellType, SubCellType, CellType_predict, GSE, Summary, Citation, PMID, cluster) %>%
#                        summarise(cell_ct = n())
#dbWriteTable(pool, 'meta_only_grouped_stats', grouped_stats, overwrite = TRUE)

args[8] <- Sys.Date()
dbWriteTable(pool, 'input_data', args %>% enframe(), overwrite = TRUE)
poolClose(pool)
