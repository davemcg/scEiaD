library(pool)
library(RSQLite)
library(Seurat)
library(tidyverse)
library(dtplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

git_dir = Sys.getenv('SCIAD_GIT_DIR')

load(args[1])
load(args[2])
load(args[3])
load(args[4])

# add study meta
study_meta <- read_tsv(paste0(git_dir, '/data/GEO_Study_Level_Metadata.tsv'))
umap <- umap %>% left_join(phate_2D$embedding %>% as_tibble(rownames = 'Barcode'))
umap <- umap %>% left_join(., study_meta, by = c('study_accession'))

# make metadata
metadata <- umap %>% select(Barcode, UMAP_1, UMAP_2, PHATE1, PHATE2, nCount_RNA, nFeature_RNA, Phase, percent_mt = `percent.mt`, S_Score = `S.Score`, G2M_Score = `G2M.Score`, batch, sample_accession, study_accession, Age, library_layout, retina_region, strain, sex, Organ, Tissue, organism, Platform, UMI, Covariate, CellType, CellType_predict, TabulaMurisCellType, TabulaMurisCellType_predict, SubCellType, GSE, Summary, Citation, PMID, Design) %>% 
		mutate(SubCellType = gsub('p_','', SubCellType)) %>% mutate(SubCellType = gsub('f_','', SubCellType)) %>%
			left_join(., meta, by = 'Barcode')

metadata$SubCellType[grepl('^RB|Rods|Peri|^MG$|^Mic$', metadata$SubCellType)] <- NA
metadata$CellType[grepl('Doub|\\/Margin\\/Periocular', metadata$CellType)] <- NA


colnames(metadata)[ncol(metadata)] <- 'subcluster'
colnames(metadata)[(ncol(metadata) - 1)] <- 'cluster'

metadata <- metadata %>%  
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

# cpm <- RelativeCounts(integrated_obj@assays$RNA@counts[, umap$Barcode], scale.factor= 1e6)
# dropping cpm for raw counts for now
# cpm was causing too many viz oddities
counts <- integrated_obj@assays$RNA@counts[, umap$Barcode]

chunk_num = 20

print(dim(counts))	
# replace with new gene names including HGNC
#gene_id_converter <- read_tsv('references/ensembl_biomart_human2mouse_macaque.tsv', skip = 1,
#                              col_names= c('hs_gene_id','hs_gene_id_v', 'mm_gene_id', 'mf_gene_id',
#                                           'hs_gene_name', 'mf_gene_name', 'mm_gene_name')) %>%
#  select(-hs_gene_id_v)
gene_id_converter <- read_tsv(glue::glue('{git_dir}/data/ensembl_biomart_human2mouse_macaque_chick_ZF.tsv.gz'), skip = 1,
                              col_names= c('hs_gene_id','hs_gene_id_v',
                                           'gg_gene_id', 'gg_gene_name',
                                           'mf_gene_id', 'mf_gene_name',
                                           'mm_gene_id', 'mm_gene_name',
                                           'zf_gene_id', 'zf_gene_name',
                                           'hs_gene_name')) %>%
  select(-hs_gene_id_v)

gene_id_converter2 <- bind_rows(gene_id_converter %>% filter(!is.na(hs_gene_id)), gene_id_converter %>% filter(!is.na(mm_gene_id)) %>% mutate(hs_gene_id = mm_gene_id, hs_gene_name = mm_gene_name))
# fill in missing mouse gene id
hs_gtf <- rtracklayer::import('references/gtf/hs-homo_sapiens_anno.gtf.gz') %>% as_tibble() %>% select(hs_gene_id = gene_id, hs_gene_name = gene_name) %>% mutate(hs_gene_id = gsub('\\.\\d+','', hs_gene_id)) %>% unique()
mf_gtf <- rtracklayer::import('references/gtf/mf-macaca_mulatta_anno.gtf.gz') %>% as_tibble() %>% select(hs_gene_id = gene_id, hs_gene_name = gene_name) %>% mutate(hs_gene_id = gsub('\\.\\d+','', hs_gene_id)) %>% unique()
mm_gtf <- rtracklayer::import('references/gtf/mm-mus_musculus_anno.gtf.gz') %>% as_tibble() %>% select(hs_gene_id = gene_id, hs_gene_name = gene_name) %>% mutate(hs_gene_id = gsub('\\.\\d+','', hs_gene_id)) %>% unique()
gene_id_converter2 <- bind_rows(hs_gtf, mf_gtf, mm_gtf) %>% unique()


row.names(counts) <- row.names(counts) %>% enframe(value = 'hs_gene_id') %>% dplyr::select(-name) %>% left_join(gene_id_converter2 %>% select(hs_gene_id, hs_gene_name) %>% unique()) %>% mutate(nname = paste0(hs_gene_name, ' (', hs_gene_id, ')')) %>% pull(nname)

# for scEiaD
genes <- row.names(counts) %>% enframe(value = 'Gene') %>% dplyr::select(-name) %>% arrange()

long_data <- list()

if (args[6] == 'TRUE'){
	seurat_meta <- integrated_obj@meta.data
	seurat_meta <- seurat_meta %>% as_tibble(rownames = 'Barcode') %>% left_join(., umap %>% select(Barcode, organism, CellType_predict), by = 'Barcode')
	for (i in seq(1:chunk_num)){
		chunk <- round(nrow(counts)/chunk_num)
   		end = chunk * i
    	start = (chunk * i) - chunk + 1
    	if (i == chunk_num){
        	end = nrow(counts)
    	}

    	print(paste(start, end))
		design = model.matrix(~ as.factor(seurat_meta$batch) + as.factor(seurat_meta$CellType_predict))
								#as.numeric(seurat_meta$`percent.mt`) + 
								#as.numeric(seurat_meta$`nCount_RNA`) + 
								#as.factor(seurat_meta$organism) +
								#as.factor(seurat_meta$CellType_predict))
    	cor_mat <- limma::removeBatchEffect(VISION:::matLog2(counts[start:end,]), 
							#batch = as.factor(seurat_meta$batch),
							design = design)
		long_data[[i]] <- cor_mat %>%
    		as_tibble(rownames = 'Gene') %>%
    		pivot_longer(cols = 2:(ncol(cor_mat) + 1), names_to = 'Barcode', values_to = 'counts') %>%
    		filter(counts > 2.5)
		rm(cor_mat)
		}
	rm(seurat_meta)
} else {
	# cpm <- VISION:::matLog2(cpm)
	for (i in seq(1:chunk_num)){
		chunk <- round(ncol(counts)/chunk_num)
		end = chunk * i
 		start = (chunk * i) - chunk + 1
 		if (i == chunk_num){
			end = ncol(counts)
		}
		
		print(paste(start, end))
		long_data[[i]] <- counts[,start:end] %>%
			as.matrix() %>%
			as_tibble(rownames = 'Gene') %>%
			pivot_longer(cols = 2:(end - start + 2), names_to = 'Barcode', values_to = 'counts') %>%
			filter(counts > 0) %>% 
			mutate(counts = log2(counts + 1))
	}
}
rm(integrated_obj)

long <- bind_rows(long_data)
pool <- dbPool(RSQLite::SQLite(), dbname = args[5])
dbWriteTable(pool, "counts", long, overwrite = TRUE)
db_create_index(pool, table = 'counts', columns = c('Gene'))
db_create_index(pool, table = 'counts', columns = c('Barcode'))




dbWriteTable(pool, "metadata", metadata, overwrite = TRUE)
db_create_index(pool, table = 'metadata', columns = c('Barcode'))

dbWriteTable(pool, "metadata_filter", meta_filter, overwrite = TRUE)
db_create_index(pool, table = 'metadata_filter', columns = c('Barcode'))

dbWriteTable(pool, "genes", genes, overwrite = TRUE)

# dtplyr
long <- lazy_dt(long)
meta_filter <- lazy_dt(meta_filter)

long <- left_join(long, meta_filter %>% select(Barcode, batch, study_accession, library_layout, retina_region, strain, sex, Organ, Tissue, organism,Stage, Platform, Covariate, CellType, SubCellType, CellType_predict, TabulaMurisCellType, TabulaMurisCellType_predict, Phase, GSE, Citation, PMID, cluster), by = 'Barcode') 

grouped_stats <- long %>% group_by(Gene, batch, study_accession, library_layout, retina_region, strain, sex, Organ, Tissue, organism,Stage, Platform, Covariate, CellType, SubCellType, CellType_predict, TabulaMurisCellType, TabulaMurisCellType_predict, Phase, GSE, Citation, PMID, cluster) %>%
						summarise(cell_ct = n(), cell_exp_ct = sum(counts > 0), counts = mean(counts)) 
# turn back into tibble
# previous left_join and group_by were lazily evaluated
grouped_stats <- as_tibble(grouped_stats)

dbWriteTable(pool, 'grouped_stats', grouped_stats, overwrite = TRUE)
db_create_index(pool, table = 'grouped_stats', columns = c('Gene'))

#meta_only_grouped_stats <- long %>% group_by(batch, study_accession, library_layout, retina_region, strain, sex, Organ, Tissue, organism,Stage, Platform, Covariate, CellType, SubCellType, CellType_predict, GSE, Summary, Citation, PMID, cluster) %>%
#                        summarise(cell_ct = n())
#dbWriteTable(pool, 'meta_only_grouped_stats', grouped_stats, overwrite = TRUE)

args[7] <- Sys.Date()
dbWriteTable(pool, 'input_data', args %>% enframe(), overwrite = TRUE)
poolClose(pool)
