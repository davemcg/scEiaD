library(tidyverse)
library(Matrix)
args = commandArgs(trailingOnly=TRUE)

organism <- args[1]
if (organism == 'Macaca_fascicularis') {
	load('quant/Macaca_fascicularis/hs_full_sparse_matrix.Rdata') # human fasta based counts
	hs_fasta <- m
	load('quant/Macaca_fascicularis/mf_full_sparse_matrix.Rdata') # macaque fasta based counts
	mf_fasta <- m
	
	# calc sum expression of each gene
	## for human fasta
	hs_rowSums <- rowSums(hs_fasta)
	names(hs_rowSums) <- row.names(hs_fasta)
	## for macaque fasta
	mf_rowSums <- rowSums(mf_fasta)
	names(mf_rowSums) <- row.names(mf_fasta)
	
	
	joined <- left_join(hs_rowSums %>% enframe() %>% rename(hs_val = value), 
						mf_rowSums %>% enframe() %>% rename(mf_val = value), by = 'name')
	joined[is.na(joined)] <- 0
	
	hs_genes <- joined %>% filter(hs_val > mf_val) %>% pull(name)
	mf_genes <- joined %>% filter(hs_val <= mf_val) %>% pull(name)
	
	# awkward two step where we take genes that ARE NOT in macaque and put back in the human list
	hs_genes <- c(hs_genes, mf_genes[!mf_genes %in% row.names(mf_fasta)])
	mf_genes <- mf_genes[mf_genes %in% row.names(mf_fasta)]
	
	# the human-based and macaque-based matrixces have different cells
	# this is because empty cells get dropped 
	mf_specific_cells <- colnames(mf_fasta)[!colnames(mf_fasta) %in% colnames(hs_fasta)]
	hs_specific_cells <- colnames(hs_fasta)[!colnames(hs_fasta) %in% colnames(mf_fasta)]
	
	# build empty sparse matrices
	mf_specific_matrix <- Matrix(0, nrow = nrow(hs_fasta), ncol = length(mf_specific_cells), sparse = TRUE)
	colnames(mf_specific_matrix) <- mf_specific_cells
	row.names(mf_specific_matrix) <- row.names(hs_fasta)
	
	hs_specific_matrix <- Matrix(0, nrow = nrow(mf_fasta), ncol = length(hs_specific_cells), sparse = TRUE)
	colnames(hs_specific_matrix) <- hs_specific_cells
	row.names(hs_specific_matrix) <- row.names(mf_fasta)
	
	# append missing cell (with 0 vals) to existing matrices
	hs_fasta1 <- cbind(hs_fasta, mf_specific_matrix)
	mf_fasta1 <- cbind(mf_fasta, hs_specific_matrix)
	
	# reorder columns to match up hs with mf
	hs_fasta1 <- hs_fasta1[,colnames(mf_fasta1)]
	
	new_mat <- rbind(hs_fasta1[hs_genes,],
					 mf_fasta1[mf_genes,]) 
	
	# sort rows by gene name
	new_mat <- new_mat[row.names(new_mat) %>% sort(),]
	
	# remove empty rows
	m <- new_mat[rowSums(new_mat) > 0, ]


	# read in cell info and combine hs and mf
	mf_cell <- read_tsv('Macaca_fascicularis_mf_cell_info.tsv')
	hs_cell <- read_tsv('Macaca_fascicularis_hs_cell_info.tsv')

	new_cell_info <- bind_rows(mf_cell %>% select(-name), hs_cell %>% select(-name)) %>% unique()
	
}	else {
	print(paste0('Copying ', organism))
	load(args[5])
	new_cell_info <- read_tsv(args[4]) %>% select(-name)
}

save(m, file = args[2], compress = FALSE)
write_tsv(new_cell_info, path = args[3])
