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

def SRS_info(SRS, data_to_return):
	organism = SRS_dict[SRS]['organism']
	if organism.lower() == 'mus_musculus':
		idx = 'references/kallisto_idx/gencode.vM22.pc_transcripts.fa.gz.idx'
		txnames = 'references/gencode.vM22.metadata.MGI_tx_mapping.tsv'
	elif organism.lower() == 'homo_sapiens':
		idx = 'references/kallisto_idx/gencode.v31.pc_transcripts.fa.gz.idx'
		txnames = 'references/gencode.v31.metadata.HGNC_tx_mapping.tsv'
	elif organism.lower() == 'macaca_fascicularis':
		idx = 'references/kallisto_idx/GCF_000364345.1_Macaca_fascicularis_5.0_rna.fna.gz.idx'
		txnames = 'references/GCF_000364345.1_Macaca_fascicularis_5.0_tx_mapping.tsv'
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

method = ['CCA', 'scanorama', 'harmony', 'fastMNN', 'combat', 'none']
transform = ['standard', 'SCT']
covariate = ['study_accession', 'batch']
organism = ['Mus_musculus', 'Macaca_fascicularis', 'Homo_sapiens']

wildcard_constraints:
	SRS = '|'.join(SRS_UMI_samples + SRS_nonUMI_samples),
	method = '|'.join(method),
	transform = '|'.join(transform),
	covariate = '|'.join(covariate),
	organism = '|'.join(organism)

rule all:
	input:
		# study_accession - early has only one covariate, so skip
		expand('plots/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.color_study__facet_age.pdf', \
				transform = transform, \
				method = method, \
				organism = 'Mus_musculus', \
				partition = ['late', 'full'], \
				covariate = covariate, \
				dims = [20,50]),
		expand('plots/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.color_study__facet_age.pdf', \
				transform = transform, \
				method = method, \
				organism = 'Mus_musculus', \
				partition = ['early'], \
				covariate = ['batch'], \
				dims = [20,50]),
		expand('plots/well_supported_celltypes/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.WellSupportedCells.color_study__facet_celltype.pdf', \
				transform = transform, \
				method = method, \
				organism = 'Mus_musculus', \
				partition = ['full'], \
				covariate = covariate, \
				dims = [20,50])
		#'quant/Mus_musculus/scTransformCCA_merged_Embryonic.seuratV3.Rdata',
		#'quant/Mus_musculus/scTransformCCA_merged_Postnatal.seuratV3.Rdata',
		#expand('quant/{SRS}/genecount/matrix.Rdata', SRS = SRS_UMI_samples), # UMI data
		#expand('quant/{SRS}/abundance.tsv.gz', SRS = SRS_nonUMI_samples) # non UMI data
		# expand('quant/{SRS}/output.bus', SRS = SRS_UMI_samples)

# mouse, human, macaque fasta and gtf
rule download_references:
	output:
		mouse_fasta = 'references/gencode.vM22.pc_transcripts.fa.gz',
		macaque_fasta = 'references/GCF_000364345.1_Macaca_fascicularis_5.0_rna.fna.gz',
		human_fasta = 'references/gencode.v31.pc_transcripts.fa.gz'
	shell:
		"""
		mkdir -p references
		wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.pc_transcripts.fa.gz  
		wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_rna.fna.gz
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
		mf = 'references/GCF_000364345.1_Macaca_fascicularis_5.0_tx_mapping.tsv',
		hs = 'references/gencode.v31.metadata.HGNC_tx_mapping.tsv',
		mm = 'references/gencode.vM22.metadata.MGI_tx_mapping.tsv'
	shell:
		"""
		zgrep "^>" references/gencode.vM22.pc_transcripts.fa.gz | \
			sed 's/>//g' | \
			awk 'BEGIN {{OFS = "\t"; FS = "|"}}; {{print $0, $2, $6}}' > {output.mm}

		wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_feature_table.txt.gz
		zcat GCF_000364345.1_Macaca_fascicularis_5.0_feature_table.txt.gz | \
			grep "^mRNA" | \
			awk -v OFS="\t" 'BEGIN {{FS="\t"}}; {{print $11,$15,$14}}' > {output.mf}
		rm GCF_000364345.1_Macaca_fascicularis_5.0_feature_table.txt.gz

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


rule make_seurat_objs:
	input:
		srr_metadata = config['srr_sample_file'],		
		tx_map = lambda wildcards: SRS_info(organism_well_dict[wildcards.organism][0], 'tx'),
		well = 'quant/{organism}/counts.Rdata',
		droplet = lambda wildcards: expand('quant/{SRS}/genecount/matrix.Rdata', SRS = organism_droplet_dict[wildcards.organism])
	output:
		#cell_info = '{organism}_cell_info.Rdata',
		seurat = 'seurat_obj/{organism}__{transform}__{partition}__{covariate}.seuratV3.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/build_seurat_obj_classic.R {output.seurat} {wildcards.partition} {wildcards.covariate} {wildcards.transform} {input}	
		"""



rule integrate:
	input:
		'seurat_obj/{organism}__{transform}__{partition}__{covariate}.seuratV3.Rdata'
	output:
		('seurat_obj/{organism}__{transform}__{partition}__{covariate}__{method}.seuratV3.Rdata')
	run:
		if wildcards.method != 'scanorama':
			job = "module load R/3.6; \
					Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/merge_methods.R \
					  {method} {transform} {covariate} {input} {output}".format(method = wildcards.method, \
																	transform = wildcards.transform, \
																    covariate = wildcards.covariate, \
											                        output = output, input = input)
		else:
			job = "source /data/mcgaugheyd/conda/etc/profile.d/conda.sh; \
					conda init bash; conda activate scanorama; \
					Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/merge_methods.R \
					  {method} {transform} {covariate} {input} {output}".format(method = wildcards.method, \
																	transform = wildcards.transform, \
																	covariate = wildcards.covariate, \
																	output = output, input = input)
		sp.run("echo " + job + '\n', shell = True)
		sp.run(job, shell = True)

#rule label_known_cells_with_type:
#	input:
#		'seurat_obj/{organism}__standard_and_SCT__{partition}__{covariate}.seuratV3.Rdata'
#	output:
#		'{organism}_cell_info_labelled.Rdata'
#	shell:
#		"""
#		module load R/3.6
#		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/label_known_cells.R {wildcards.organism}_cell_info.Rdata {output}
#		"""

rule calculate_umap_and_cluster:
	input:
		'{organism}_cell_info_labelled.Rdata',
		'seurat_obj/{organism}__{transform}__{partition}__{covariate}__{method}.seuratV3.Rdata'
	output:
		'seurat_obj/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.umap.seuratV3.Rdata',
		'umap/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.umap.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/calculate_umap_and_cluster.R {wildcards.method} {wildcards.dims} {input} {output}
		"""

rule plot_integration:
	input:
		'umap/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.umap.Rdata'
	output:
		'plots/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.color_study__facet_age.pdf',
		'plots/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.color_study__facet_batch.pdf',
		'plots/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.color_paper__facet_celltype.pdf',
		'plots/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.color_celltype__facet_cluster.pdf'
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/plot_assess_integration.R {input} {output}
		"""

rule plot_integration_with_well_supported_cell_types:
	input:
		'umap/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.umap.Rdata'
	output:
		'plots/well_supported_celltypes/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.WellSupportedCells.color_study__facet_celltype.pdf',
		'plots/well_supported_celltypes/{organism}__{transform}__{partition}__{covariate}__{method}__dims{dims}.WellSupportedCells.color_celltype.pdf',
	shell:
		"""
		module load R/3.6
		Rscript /home/mcgaugheyd/git/massive_integrated_eye_scRNA/src/plot_assess_integration_specific_types.R {input} {output}
		"""
