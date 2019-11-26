import subprocess as sp

srr_sample_file = config['srr_sample_file']
#srr_sample_discrepancy_file = config['srr_discrepancy_file']

# builds dictionary of dictionaries where first dict key is SRS 
# and second dict key are SRS properties
def metadata_builder(file, SRS_dict = {}, discrepancy = False):
	with open(file) as file:
		for line in file:
			info = line.strip('\n').split('\t')
			if info[0] == 'sample_accession':
				continue
			SRS = info[0]
			if SRS not in SRS_dict:
				SRS_dict[SRS]={'SRR': [info[1]],
					    	  'paired':True if info[2]=='PAIRED' else False, 
					          'organism':info[3].replace(' ', '_'),
		            	      'tech':info[4],
						      'UMI':True if info[5]=='YES' else False,
							  'Study': info[6]}
			else:
				# this is mostly for SRA having the 'paired' status wrong
				# don't want to hand-edit the main metadata file
				# so I think better to create a new file with
				# hand edited values for just the ones i want to change
				if discrepancy:
					runs = SRS_dict[SRS]['SRR']
					SRS_dict[SRS] = {'SRR':runs,
									 'paired':True if info[2]=='PAIRED' else False,
									 'organism':info[3],
									 'tech':info[4],
									 'UMI':True if info[5]=='YES' else False,
									 'Study': info[6]}
				else:
					runs = SRS_dict[SRS]['SRR']
					runs.append(info[1])
					SRS_dict[SRS]['SRR'] = runs
	return(SRS_dict)

SRS_dict = metadata_builder(srr_sample_file)
# hand edited file which corrects mistakes that in the 
# various databases
# SRS_dict = metadata_builder(srr_sample_discrepancy_file, SRS_dict, discrepancy = True)

# build organism <-> SRS dict for nonUMI data
organism_well_dict = {}
for x in SRS_dict:
	if not SRS_dict[x]['UMI'] or not SRS_dict[x]['paired']:
		organism = SRS_dict[x]['organism']
		if organism not in organism_well_dict:
			organism_well_dict[organism] = [x]
		else:
			srs = organism_well_dict[organism]
			srs.append(x)
			organism_well_dict[organism] = srs
# build organsim <-> SRS dict for UMI/droplet data
organism_droplet_dict = {}
for x in SRS_dict:
	if SRS_dict[x]['UMI'] and SRS_dict[x]['paired']:
		organism = SRS_dict[x]['organism']
		if organism not in organism_droplet_dict:
			organism_droplet_dict[organism] = [x]
		else:
			srs = organism_droplet_dict[organism]
			srs.append(x)
			organism_droplet_dict[organism] = srs

def lookup_run_from_SRS(SRS):
	SRR_files=SRS_dict[SRS]['SRR']
	out = []
	for SRR in SRR_files:
		if SRS_dict[SRS]['paired']:
			#PE
			out.append('fastq/{}_1.fastq.gz'.format(SRR))
			out.append('fastq/{}_2.fastq.gz'.format(SRR))
		else:
			#SE
			out.append('fastq/{}.fastq.gz'.format(SRR))
	return(out)

# return dummy well file for macaca ,as there is no well data at this time
def well_and_droplet_input(organism):
	if organism == 'Macaca_fascicularis':
		out = ['quant/' + x + '/genecount/matrix.Rdata' for x in organism_droplet_dict[organism]]
	else:
		out = ['quant/' + organism + '/counts.Rdata']	+ \
				['quant/' + x + '/genecount/matrix.Rdata' for x in organism_droplet_dict[organism]]
	return(out)

def SRS_info(SRS, data_to_return):
	organism = SRS_dict[SRS]['organism']
	if organism.lower() == 'mus_musculus':
		idx = 'references/kallisto_idx/gencode.vM22.pc_transcripts.fa.gz.idx'
		txnames = 'references/gencode.vM22.metadata.MGI_tx_mapping.tsv'
	elif organism.lower() == 'homo_sapiens':
		idx = 'references/kallisto_idx/gencode.v31.pc_transcripts.fa.gz.idx'
		txnames = 'references/gencode.v31.metadata.HGNC_tx_mapping.tsv'
	elif organism.lower() == 'macaca_fascicularis':
		idx = 'references/kallisto_idx/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz.idx'
		txnames = 'references/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all_tx_mapping.tsv'
	else:
		print(SRS + ' ' + organism + " NO SPECIES MATCH!")
	if data_to_return == 'idx':
		out = idx
	else:
		out = txnames
	return(out)


SRS_UMI_samples = []
SRS_nonUMI_samples = []
for SRS in SRS_dict.keys():
	if SRS_dict[SRS]['UMI'] and SRS_dict[SRS]['paired']:
		SRS_UMI_samples.append(SRS)
	elif SRS_dict[SRS]['tech'] != 'BULK':
		SRS_nonUMI_samples.append(SRS)

method = ['scVI','CCA', 'scanorama', 'harmony', 'fastMNN', 'combat', 'none']
transform = ['counts','standard', 'SCT','scran']
covariate = ['study_accession', 'batch']
organism = ['Mus_musculus', 'Macaca_fascicularis', 'Homo_sapiens']
combination = ['Mus_musculus', 'Mus_musculus_Macaca_fascicularis', 'Mus_musculus_Macaca_fascicularis_Homo_sapiens']
dims = [8,10,20,25,30,50,75,100,200]

wildcard_constraints:
	SRS = '|'.join(SRS_UMI_samples + SRS_nonUMI_samples),
	method = '|'.join(method),
	transform = '|'.join(transform),
	covariate = '|'.join(covariate),
	organism = '|'.join(organism),
	nfeatures = '|'.join([str(x) for x in [2000,5000]]),
	dims = '|'.join([str(x) for x in dims])
	

rule all:
	input:
		'references/kallisto_idx/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz.idx',
		'quant/Macaca_fascicularis/full_sparse_matrix.Rdata',
		'seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__scran__full__batch.seuratV3.Rdata',
		'seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__SCT__full__batch.seuratV3.Rdata',
		'seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__standard__full__batch.seuratV3.Rdata',
		expand('quant/{organism}/full_sparse_matrix.Rdata', organism = organism),
		#expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.big_plot.png', \
		#		transform = ['scran', 'standard'], \
		#		method = ['fastMNN'], \
		#		combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
		#		partition = ['full'], \
		#		n_features = [2000], \
		#		covariate = ['batch'], \
		#		dims = [30,50,100,200],
		#		dist = [0.001,0,1, 0.3],
		#		neighbors = [5, 30, 50]),
		expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.big_plot.png', \
				transform = ['counts'], \
				method = ['scVI'], \
				combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
				partition = ['full'], \
				n_features = [2000, 5000], \
				covariate = ['batch'], \
				dims = [8,10,20,30,50,100,200],
				dist = [0.001,0.1, 0.3],
				neighbors = [5, 15, 30, 50, 100]),
		#expand('cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn{knn}.cluster.Rdata', \
		#		transform = ['scran', 'standard'], \
		#		method = ['fastMNN'], \
		#		combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
		#		partition = ['full'], \
		#		n_features = [2000], \
		#		covariate = ['batch'], \
		#		knn = [5, 7, 10], \
		#		dims = [30,50,100,200]),
		expand('cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn{knn}.cluster.Rdata', \
				transform = ['counts'], \
				method = ['scVI'], \
				n_features = [2000,5000], \
				combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
				partition = ['full'], \
				covariate = ['batch'], \
				knn = [4, 5, 7, 10], \
				dims = [8,10,20,30,50,100,200]),
		expand('quant/{SRS}/abundance.tsv.gz', SRS = SRS_nonUMI_samples), # non UMI data
		expand('quant/{SRS}/output.bus', SRS = SRS_UMI_samples)

# mouse, human, macaque fasta and gtf
rule download_references:
	output:
		mouse_fasta = 'references/gencode.vM22.pc_transcripts.fa.gz',
		macaque_fasta = 'references/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz',
		human_fasta = 'references/gencode.v31.pc_transcripts.fa.gz'
	shell:
		"""
		mkdir -p references
		wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.pc_transcripts.fa.gz  
		#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_rna.fna.gz
		wget ftp://ftp.ensembl.org/pub/release-98/fasta/macaca_fascicularis/cdna/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz
		wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.pc_transcripts.fa.gz
		mv *fa*gz references/
		"""

# need to make mouse, human, macaque
rule kallisto_index:
	input:
		'references/{fasta}'
	output:
		'references/kallisto_idx/{fasta}.idx'
	shell:
		"""
		kallisto index {input} -i {output}
		"""

# get / make the bustool count tx file
# tsv three column
# transcript_ID gene_ID gene_name
# ENST ENSG gene_name
rule tx_gene_mapping:
	input:
		'references/gencode.vM22.pc_transcripts.fa.gz'
	output:
		mf = 'references/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all_tx_mapping.tsv',
		hs = 'references/gencode.v31.metadata.HGNC_tx_mapping.tsv',
		mm = 'references/gencode.vM22.metadata.MGI_tx_mapping.tsv'
	shell:
		"""
		zgrep "^>" references/gencode.vM22.pc_transcripts.fa.gz | \
			sed 's/>//g' | \
			awk 'BEGIN {{OFS = "\t"; FS = "|"}}; {{print $0, $2, $6}}' > {output.mm}

		zgrep "^>" references/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz | \
			sed 's/>//g' > mf.header
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/macaque_ensembl_fasta_header.R mf.header {output.mf}
		rm mf.header
		#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_feature_table.txt.gz
		#zcat GCF_000364345.1_Macaca_fascicularis_5.0_feature_table.txt.gz | \
		#	grep "^mRNA" | \
		#	awk -v OFS="\t" 'BEGIN {{FS="\t"}}; {{print $11,$15,$14}}' > {output.mf}
		#rm GCF_000364345.1_Macaca_fascicularis_5.0_feature_table.txt.gz

		zgrep "^>" references/gencode.v31.pc_transcripts.fa.gz | \
			sed 's/>//g' | \
			awk 'BEGIN {{OFS = "\t"; FS = "|"}}; {{print $0, $2, $6}}' > {output.hs}
		"""

# this does the pseudoalignment for UMI data (e.g. 10x)
rule kallisto_bus:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS),
		idx = lambda wildcards: SRS_info(wildcards.SRS, 'idx')
	output:
		bus = 'quant/{SRS}/output.bus',
		ec = 'quant/{SRS}/matrix.ec',
		tx_name = 'quant/{SRS}/transcripts.txt'
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech'],
		paired = lambda wildcards: SRS_dict[wildcards.SRS]['paired']
	run:
		if params.paired:
			job = "kallisto bus -x {tech} \
					-i {idx} -o quant/{SRS} {fastq}".format(fastq = input.fastq,
                                                    tech = params.tech,
													idx = input.idx,
													SRS = wildcards.SRS)
		else:
			job = "kallisto bus --single -x {tech} \
					-i {idx} -o quant/{SRS} {fastq}".format(fastq = input.fastq,
                                                    tech = params.tech,
													idx = input.idx,
													SRS = wildcards.SRS)
		sp.run("echo " + job + '\n', shell = True)
		sp.run(job, shell = True)	

# pseudoaligment for nonUMI data (e.g. smartseq)
rule kallisto_quant:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS),
		idx = lambda wildcards: SRS_info(wildcards.SRS, 'idx')
	output:
		quant = 'quant/{SRS}/abundance.tsv.gz'
	params:
		paired = lambda wildcards: SRS_dict[wildcards.SRS]['paired']
	threads: 8
	run:
		if params.paired:
			job = "kallisto quant -t {t} -b 100 --plaintext --bias \
					-i {idx} -o quant/{SRS} {fastq}".format(fastq = input.fastq,
                                                    t = threads,
													idx = input.idx,
													SRS = wildcards.SRS)
		else:
			job = "kallisto quant --single -l 200 -s 30  -t {t} -b 100 --plaintext --bias \
					-i {idx} -o quant/{SRS} {fastq}".format(fastq = input.fastq,
                                                    t = threads,
													idx = input.idx,
													SRS = wildcards.SRS)
		sp.run("echo " + job + '\n', shell = True)
		sp.run(job, shell = True)	
		sp.run('gzip quant/' + wildcards.SRS + '/*tsv', shell = True)
			
			
# sorting required for whitelist creation and correction
# make these temp files
rule bustools_sort:
	input:
		'quant/{SRS}/output.bus'
	output:
		('quant/{SRS}/output.sorted.bus')
	threads: 4
	shell:
		"""
		/home/mcgaugheyd/git/bustools/build/src/./bustools sort -t {threads} -m 16G \
			{input} \
			-o {output}
		"""

# find barcodes, correct barcodes
# make these temp files
rule bustools_whitelist_correct_count:
	input:
		bus = 'quant/{SRS}/output.sorted.bus',
		ec = 'quant/{SRS}/matrix.ec',
		tx_name = 'quant/{SRS}/transcripts.txt',
		tx_map = lambda wildcards: SRS_info(wildcards.SRS, 'tx')
	output:
		whitelist = 'whitelist/{SRS}_whitelist',
		bus_matrix = 'quant/{SRS}/genecount/gene.mtx'
	params:
		bus_out = 'quant/{SRS}/genecount/gene'
	shell:
		"""
		/home/mcgaugheyd/git/bustools/build/src/./bustools whitelist \
			{input.bus} \
			-o {output.whitelist}

		/home/mcgaugheyd/git/bustools/build/src/./bustools correct \
			{input.bus} \
			-w {output.whitelist} \
			-p | \
		/home/mcgaugheyd/git/bustools/build/src/./bustools count \
			-o {params.bus_out} \
			-g {input.tx_map} \
			-e {input.ec} \
			-t {input.tx_name} \
			--genecounts -
		"""

rule create_sparse_matrix:
	input:
		'quant/{SRS}/genecount/gene.mtx'
	output:
		stats = 'quant/{SRS}/genecount/stats.tsv',
		matrix = 'quant/{SRS}/genecount/matrix.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/remove_empty_UMI_make_sparse_matrix.R {wildcards.SRS} {output.matrix} {output.stats}
		"""		

rule merge_nonUMI_quant_by_organism:
	input:
		quant = lambda wildcards: expand('quant/{SRS}/abundance.tsv.gz', SRS = organism_well_dict[wildcards.organism]),
		tx_map = lambda wildcards: SRS_info(organism_well_dict[wildcards.organism][0], 'tx')
	output:
		'quant/{organism}/counts.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/merge_nonUMI_quant_by_organism.R {output} {input.tx_map} {input.quant} 
		"""

rule combine_well_and_umi:
	input:
		srr_metadata = config['srr_sample_file'],
		tx_map = lambda wildcards: SRS_info(organism_droplet_dict[wildcards.organism][0], 'tx'),
		counts = lambda wildcards: well_and_droplet_input(wildcards.organism)
	output:
		cell_info = '{organism}_cell_info.tsv',
		matrix = 'quant/{organism}/full_sparse_matrix.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/build_sparse_matrix.R {wildcards.organism} {output} {input}
		"""

localrules: cat_cell_info
rule cat_cell_info:
	input:
		expand('{organism}_cell_info.tsv', \
			organism = organism)
	output:
		'cell_info.tsv'
	shell:
		"""
		cat {input} | head -n 1 > header
		cat header <( grep -v "^name" {input}) > {output}
		"""
		
rule make_seurat_objs:
	input:
		'cell_info.tsv',
		expand('quant/{organism}/full_sparse_matrix.Rdata', \
			 organism = organism)
	output:
		seurat = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}.seuratV3.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/build_seurat_obj_classic.R \
			{output.seurat} {wildcards.partition} {wildcards.covariate} {wildcards.transform} \
			{wildcards.combination} {wildcards.n_features} {input}	
		"""

rule integrate:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}.seuratV3.Rdata',
	output:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.seuratV3.Rdata'
	run:
		job = "module load R/3.6; \
				Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/merge_methods.R \
				  {method} {transform} {covariate} {dims} {input} {output}".format(\
						method = wildcards.method, \
						transform = wildcards.transform, \
					    covariate = wildcards.covariate, \
						dims = wildcards.dims, \
			            output = output, input = input)
		sp.run("echo \"" +  job + "\"\n", shell = True)
		sp.run(job, shell = True)


rule label_known_cells_with_type:
	input:
		'cell_info.tsv',
		config['srr_sample_file']		
	output:
		'cell_info_labelled.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/label_known_cells.R 
		"""

rule predict_missing_cell_types:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.seuratV3.Rdata',
		'cell_info_labelled.Rdata'
	output:
		'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}_cell_info_predictions.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/transfer_labels.R {input} {wildcards.transform} {output}
		"""
		
rule calculate_umap:
	input:
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.seuratV3.Rdata'
	output:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/calculate_umap_and_cluster.R \
			{wildcards.method} {wildcards.dims} {wildcards.dist} {wildcards.neighbors} 1 FALSE TRUE {input} {output}
		"""

rule calculate_cluster:
	input:
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.seuratV3.Rdata'
	output:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn{knn}.cluster.seuratV3.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/calculate_umap_and_cluster.R \
			{wildcards.method} {wildcards.dims} 1 1 {wildcards.knn} TRUE FALSE {input} {output}
		"""

rule extract_umap:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn5.cluster.Rdata',
		'cell_info_labelled.Rdata',
		'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}_cell_info_predictions.Rdata'
	output:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/extract_umap.R \
			{input} {output} {wildcards.method}	
		"""

rule extract_cluster:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn{knn}.cluster.seuratV3.Rdata'
	output:
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn{knn}.cluster.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__knn{knn}.graph.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/extract_cluster.R \
			{input} {output}
		"""




rule plot_integration:
	input:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	output:
		'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.big_plot.png'
		#'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.color_study__facet_age.pdf',
		#'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.color_study__facet_batch.pdf',
		#'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.color_paper__facet_celltype.pdf',
		#'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.color_celltype__facet_cluster.pdf'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/big_plots.R {input} {output}
		"""

rule plot_integration_with_well_supported_cell_types:
	input:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	output:
		'plots/well_supported_celltypes/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.WellSupportedCells.color_study__facet_celltype.pdf',
		'plots/well_supported_celltypes/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.WellSupportedCells.color_celltype.pdf',
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/plot_assess_integration_specific_types.R {input} {output}
		"""

rule cluster_assessment:
	input:
		expand('umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__mindist{dist}__nneighbors{neighbors}.umap.Rdata', \
				transform = transform, \
				method = method, \
				combination = 'Mus_musculus', \
				partition = ['full'], \
				n_features = [2000], \
				covariate = covariate, \
				dist = [0.001, 0.3, 0.5], \
				neighbors = [5, 30, 50],\
				dims = dims)
	output:
		'plots/well_supported_celltypes/{organism}.mean_cluster_purity_by_cell_type.pdf',
		'plots/well_supported_celltypes/{organism}.mean_cluster_purity_by_study_and_cell_type.pdf',
		'plots/well_supported_celltypes/{organism}.median_cluster_area_chull.pdf'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/cluster_purity.R umap {output}
		"""
