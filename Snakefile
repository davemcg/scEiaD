import pprint
pp = pprint.PrettyPrinter(width=41, compact=True) 
import subprocess as sp

srr_sample_file = config['srr_sample_file']
#srr_sample_discrepancy_file = config['srr_discrepancy_file']

# builds dictionary of dictionaries where first dict key is SRS 
# and second dict key are SRS properties
def metadata_builder(file, SRS_dict = {}, discrepancy = False):
	with open(file) as file:
		for line in file:
			if line[0] == '#':
				continue
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
# pp.pprint(SRS_dict)
# hand edited file which corrects mistakes that in the 
# various databases
# SRS_dict = metadata_builder(srr_sample_discrepancy_file, SRS_dict, discrepancy = True)

# build organism <-> SRS dict for nonUMI data
organism_well_dict = {}
organism_welltech_dict ={'Homo_sapiens':[], 'Mus_musculus':[], 'Macaca_fascicularis':[] }
		
	
for x in SRS_dict:
	if not SRS_dict[x]['UMI'] or not SRS_dict[x]['paired']:
		organism = SRS_dict[x]['organism']
		tech = SRS_dict[x]['tech']
		if tech not in organism_welltech_dict[organism]:
			organism_welltech_dict[organism].append(tech)
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


def lookup_run_from_SRS(SRS, fqp):
	SRR_files=SRS_dict[SRS]['SRR']
	out = []
	for SRR in SRR_files:
		if SRS_dict[SRS]['paired']:
			#PE
			out.append(f'{fqp}/fastq/{SRR}_1.fastq.gz')
			out.append(f'{fqp}/fastq/{SRR}_2.fastq.gz')
		else:
			#SE
			out.append(f'{fqp}/fastq/{SRR}.fastq.gz')
	return(out)

# return dummy well file for macaca ,as there is no well data at this time
def well_and_droplet_input(organism, reference, quant_path, SRS_dict, otd):
	if organism == 'Macaca_fascicularis':
		out = [f'{quant_path}/quant/{srs}/{SRS_dict[srs]["tech"]}/{reference}/genecount/matrix.Rdata' for srs in organism_droplet_dict[organism]]
	else:
		out = [f'{quant_path}/quant/{organism}/{tech}/{reference}__counts.Rdata' for tech in otd[organism]] + \
				[f'{quant_path}/quant/{srs}/{SRS_dict[srs]["tech"]}/{reference}/genecount/matrix.Rdata' for srs in organism_droplet_dict[organism]]
	return(out)

def REF_idx(ref, data_to_return):
	organism = SRS_dict[SRS]['organism']
	if ref == 'mm':
		idx = 'references/kallisto_idx/gencode.vM25.pc_transcripts.fa.gz.idx'
		txnames = 'references/gencode.vM25.metadata.MGI_tx_mapping.tsv'
	elif ref == 'hs-homo_sapiens':
		idx = 'references/kallisto_idx/gencode.v34.pc_transcripts.fa.gz.idx'
		txnames = 'references/gencode.v34.metadata.HGNC_tx_mapping.tsv'
	elif ref == 'mf':
		idx = 'references/kallisto_idx/Macaca_mulatta.Mmul_10.cdna.all.fa.gz.idx'
		txnames = 'references/Macaca_mulatta.Mmul_10.cdna.all_tx_mapping.tsv'
	else:
		print(ref + " NO REF MATCH!")
	if data_to_return == 'idx':
		out = idx
	else:
		out = txnames
	return(out)




def ORG_ref(organism):
	if organism.lower() == 'mus_musculus':
		out = ['mm']
	elif organism.lower() == 'homo_sapiens':
		out = ['hs']
	elif organism.lower() == 'macaca_fascicularis':
		out = ['hs','mf']
	else:
		print(organism + ' NO MATCH')
	return(out)

git_dir = config['git_dir']
bustools_path = config['bustools_path']
working_dir = config['working_dir']
conda_dir = config['conda_dir']
fastq_path = config['fastq_path']
fi_tsne_dir = config['fi_tsne_dir']
quant_path = config['quant_path']
config_abspath=config['config_abspath']

SRS_UMI_samples = []
SRS_nonUMI_samples = []
for SRS in SRS_dict.keys():
	if SRS_dict[SRS]['UMI'] and SRS_dict[SRS]['paired']:
		SRS_UMI_samples.append(SRS)
	elif SRS_dict[SRS]['tech'] != 'BULK':
		SRS_nonUMI_samples.append(SRS)

method = [ 'bbknn','insct','magic', 'scVI','CCA', 'scanorama', 'harmony', 'fastMNN', 'combat', 'none', 'liger']
transform = ['libSize', 'sqrt', 'counts','standard', 'SCT','scran']
covariate = ['study_accession', 'batch']
organism = ['Mus_musculus', 'Macaca_fascicularis', 'Homo_sapiens']
combination = ['Mus_musculus', 'Mus_musculus_Macaca_fascicularis', 'Mus_musculus_Macaca_fascicularis_Homo_sapiens', 'universe']
dims = [4,6,8,10,20,25,30,50,75,100,200]
knn = [0.2,0,4,0.6, 5, 7, 10, 15]
model = ['A', 'B', 'C', 'D', 'E', 'F', 'G'] # A is ~seuratCluster+batch+percent.mt and B is ~seuratCluster+batch+percent.mt+organism
wildcard_constraints:
	SRS = '|'.join(SRS_UMI_samples + SRS_nonUMI_samples),
	method = '|'.join(method),
	transform = '|'.join(transform),
	covariate = '|'.join(covariate),
	organism = '|'.join(organism),
	nfeatures = '|'.join([str(x) for x in [2000,5000]]),
	dims = '|'.join([str(x) for x in dims]),
	model = '|'.join(model)	

rule all:
	input:
		human_hs_quant = well_and_droplet_input('Homo_sapiens', 'hs-homo_sapiens', quant_path, SRS_dict, organism_welltech_dict)#,
		# mouse_mm_quant = well_and_droplet_input('Mus_musculus', 'mm-mus_musculus', quant_path, SRS_dict, organism_welltech_dict),
		# monkey_mf_quant = well_and_droplet_input('Macaca_fascicularis', 'mf-macaca_mulatta', quant_path, SRS_dict, organism_welltech_dict),
		# monkey_hs_quant = well_and_droplet_input('Macaca_fascicularis', 'hs-homo_sapiens', quant_path, SRS_dict, organism_welltech_dict),
		# human_dntx_quant = well_and_droplet_input('Homo_sapiens', 'DNTX', quant_path, SRS_dict, organism_welltech_dict)
		

# if config['subset_clustering'] == 'False': 
# 	rule all:
# 		input:
# 			expand(quant_path +'/quant/{sample}/hs/output.bus', sample = SRS_UMI_samples),
# 			expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plot.png', \
# 					transform = ['counts'], \
# 					method = ['scVI'], \
# 					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
# 					partition = ['universe'], \
# 					n_features = [2000, 5000, 10000], \
# 					covariate = ['batch'], \
# 					dims = [30,50,100],
# 					dist = [0.1],
# 					neighbors = [100]),
# 			expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plotPredictions.png', \
# 					transform = ['libSize'], \
# 					method = ['fastMNN'], \
# 					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
# 					partition = ['onlyWELL'], \
# 					n_features = [2000], \
# 					covariate = ['batch'], \
# 					dims = [30],
# 					dist = [0.1],
# 					neighbors = [50]),
# 			expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plot.png', \
# 					transform = ['counts'], \
# 					method = ['scVI'], \
# 					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
# 					partition = ['TabulaDroplet'], \
# 					n_features = [1000, 2000, 5000, 10000], \
# 					covariate = ['batch'], \
# 					dims = [4,6,8,10,20,30,50,100],
# 					dist = [0.001,0.1],
# 					neighbors = [15, 30, 50, 100, 500]),
# 			expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plot.png', \
# 					transform = ['sqrt','libSize','scran', 'standard', 'SCT'], \
# 					method = ['bbknn','insct',  'magic', 'scanorama', 'harmony', 'fastMNN', 'combat',  'none'], \
# 					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
# 					partition = ['TabulaDroplet'], \
# 					n_features = [2000], \
# 					covariate = ['batch'], \
# 					dims = [8, 30],
# 					dist = [0.3],
# 					neighbors = [30]),
# 			'merged_stats_2020_07_06.Rdata',
# 			#'site/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite.gz',
# else:
# 	rule all:
# 		input:	
# 			expand('plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plot.png', \
# 					transform = ['counts'], \
# 					method = ['scVI'], \
# 					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
# 					partition = ['cones', 'hc', 'rgc', 'amacrine', 'mullerglia', 'bipolar', 'rods'], \
# 					n_features = [ 1000, 2000], \
# 					covariate = ['batch'], \
# 					dims = [4,6,10,30],
# 					dist = [0.1],
# 					neighbors = [30]),
# 			'merged_stats_subsetClustering.Rdata'
	
# get annotation for each species 
rule download_annotation:
	output:
		mouse_anno='references/gtf/mm-mus_musculus_anno.gtf.gz',
		#macaque_faca_anno='references/gtf/macaca_fascicularis_anno.gtf.gz',
		macaque_mul_anno= 'references/gtf/mf-macaca_mulatta_anno.gtf.gz',
		human_anno='references/gtf/hs-homo_sapiens_anno.gtf.gz'
	shell:
		'''
		mkdir -p references
		wget -O {output.mouse_anno} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
		wget -O {output.macaque_mul_anno} ftp://ftp.ensembl.org/pub/release-99/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.99.gtf.gz
		wget -O {output.human_anno} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
		'''
		#wget -O {output.macaque_faca_anno} ftp://ftp.ensembl.org/pub/release-98/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.98.gtf.gz


rule get_velocity_files:
	input: 
		gtf = 'references/gtf/{reference}_anno.gtf.gz'
	output:
		'references/velocity/{tech}/{reference}/cDNA_introns.fa',
		'references/velocity/{tech}/{reference}/tr2g.tsv'
	params:
		out_dir = lambda wildcards: f'references/velocity/{wildcards.tech}/{wildcards.reference}/'
	shell:
		'''
		module load R/3.6
		Rscript src/get_velocity_annotation.R {input.gtf} {wildcards.reference} {wildcards.tech} {params.out_dir}
		'''



# need to make mouse, human, macaque
rule kallisto_index:
	input:
		'references/velocity/{tech}/{reference}/cDNA_introns.fa'
	output:
		'references/kallisto_idx/{tech}/{reference}.idx'
	shell:
		"""
		module load kallisto/0.46.2
		kallisto index {input} -i {output}
		"""

# get / make the bustool count tx file
# this does the pseudoalignment for UMI data (e.g. 10x)
rule kallisto_bus:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS, fastq_path),
		idx = 'references/kallisto_idx/{tech}/{reference}.idx'
	output:
		bus = quant_path + '/quant/{SRS}/{tech}/{reference}/output.bus',
		ec = quant_path + '/quant/{SRS}/{tech}/{reference}/matrix.ec',
		tx_name = quant_path + '/quant/{SRS}/{tech}/{reference}/transcripts.txt'
	threads: 4
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech'],
		paired_flag = lambda wildcards: '' if SRS_dict[wildcards.SRS]['paired'] else '--single',
		out_dir = lambda wildcards:  f'{quant_path}/quant/{wildcards.SRS}/{wildcards.tech}/{wildcards.reference}'
	shell:
		'''
		module load kallisto/0.46.2
		kallisto bus {params.paired_flag} -t 4 -x {params.tech} \
					-i {input.idx} -o {params.out_dir} {input.fastq}
		'''


# pseudoaligment for nonUMI data (e.g. smartseq)
def get_kallisto_quant_layout_flag(is_paired):
	if is_paired:
		return ''
	else: 
		return '--single -l 200 -s 30'

rule kallisto_quant:
	input:
		fastq = lambda wildcards: lookup_run_from_SRS(wildcards.SRS, fastq_path),
		idx = 'references/kallisto_idx/{tech}/{reference}.idx'
	output:
		quant = quant_path + '/quant/{SRS}/{tech}/{reference}/abundance.tsv.gz'
	params:
		paired_flag = lambda wildcards: get_kallisto_quant_layout_flag(SRS_dict[wildcards.SRS]['paired']),
		outdir =lambda wildcards:  f'{quant_path}/quant/{wildcards.SRS}/{wildcards.tech}/{wildcards.reference}'
	threads: 8
	shell:
		'''
		module load kallisto/0.46.2
		kallisto quant {params.paired_flag} -t {threads} --bias \
					-i {input.idx} -o {params.outdir} {input.fastq}
		gzip {params.outdir}/abundance.tsv
		'''
	

			
# sorting required for whitelist creation and correction
# make these temp files
rule bustools_sort:
	input:
		quant_path + '/quant/{SRS}/{tech}/{reference}/output.bus'
	output:
		temp(quant_path +'/quant/{SRS}/{tech}/{reference}/output.sorted.bus')
	threads: 4
	shell:
		"""
		{bustools_path}/./bustools sort -t {threads} -m 16G \
			{input} \
			-o {output}
		"""

# find barcodes, correct barcodes
# make these temp files
rule bustools_whitelist_correct_count:
	input:
		bus = quant_path + '/quant/{SRS}/{tech}/{reference}/output.sorted.bus',
		matrix = quant_path + '/quant/{SRS}/{tech}/{reference}/matrix.ec',
		tx_name = quant_path +'/quant/{SRS}/{tech}/{reference}/transcripts.txt',
		tx_map = 'references/velocity/{tech}/{reference}/tr2g.tsv'
	output:
		whitelist = 'whitelist/{SRS}/{tech}/{reference}_whitelist',
		spliced = quant_path +'/quant/{SRS}/{tech}/{reference}/genecount/spliced.mtx', 
		unspliced = quant_path +'/quant/{SRS}/{tech}/{reference}/genecount/unspliced.mtx', 
	params:
		bus_out =  lambda wildcards: f'{quant_path}/quant/{wildcards.SRS}/{wildcards.tech}/{wildcards.reference}/genecount/',
		vref = lambda wildcards: f'references/velocity/{wildcards.tech}/{wildcards.reference}'
	shell:
		'''
		{bustools_path}/./bustools whitelist \
			-o {output.whitelist} \
			{input.bus} 
		bustools correct -w {output.whitelist} -o {params.bus_out}/TMP.correct.sort.bus {input.bus} 
		
		bustools capture -s -x -o {params.bus_out}/TMP.spliced.bus -c {params.vref}/introns_tx_to_capture.txt -e {input.matrix} -t {params.vref}/transcripts.txt {params.bus_out}/TMP.correct.sort.bus
		bustools capture -s -x -o {params.bus_out}/TMP.unspliced.bus -c {params.vref}/cDNA_tx_to_capture.txt -e {input.matrix} -t {params.vref}/transcripts.txt {params.bus_out}/TMP.correct.sort.bus

		bustools count -o {params.bus_out}/spliced -g {params.vref}/tr2g.tsv -e {input.matrix}  -t {input.tx_name}  --genecounts {params.bus_out}/TMP.spliced.bus
		bustools count -o {params.bus_out}/unspliced -g {params.vref}/tr2g.tsv -e {input.matrix} -t {input.tx_name}  --genecounts {params.bus_out}/TMP.unspliced.bus
		rm {params.bus_out}/TMP*
		'''
	
	

## need to fix this script	
rule create_sparse_matrix:
	input:
		spliced = quant_path +'/quant/{SRS}/{tech}/{reference}/genecount/spliced.mtx', 
		unspliced = quant_path +'/quant/{SRS}/{tech}/{reference}/genecount/unspliced.mtx', 

	output:
		stats = quant_path + '/quant/{SRS}/{tech}/{reference}/genecount/stats.tsv',
		spliced_matrix = quant_path + '/quant/{SRS}/{tech}/{reference}/genecount/matrix.Rdata',
		unspliced_matrix = quant_path + '/quant/{SRS}/{tech}/{reference}/genecount/unspliced_matrix.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/remove_empty_UMI_make_sparse_matrix.R {wildcards.SRS} {wildcards.reference} {output.spliced_matrix} {output.stats} {working_dir}
		Rscript {git_dir}/src/remove_empty_UMI_make_sparse_matrix.R {wildcards.SRS} {wildcards.reference} {output.unspliced_matrix} {output.stats} {working_dir}
		"""		


rule merge_nonUMI_quant_by_organism:
	input:
		quant = lambda wildcards: expand(quant_path + '/quant/{SRS}/{{tech}}/{{reference}}/abundance.tsv.gz', 
										SRS = [srs for srs in organism_well_dict[wildcards.organism] if SRS_dict[srs]['tech'] == wildcards.tech and SRS_dict ] ),
		tx_map = 'references/velocity/{tech}/{reference}/tr2g.tsv'
	output:
		quant_path + '/quant/{organism}/{tech}/{reference}__counts.Rdata',
		quant_path + '/quant/{organism}/{tech}/{reference}__counts_tx.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/merge_nonUMI_quant_by_organism.R {output} {input.tx_map} {input.quant} 
		"""
# 144419
rule combine_well_and_umi:
	input:
		srr_metadata = config['srr_sample_file'],
		tx_map = lambda wildcards:  REF_idx(wildcards.reference, 'tx'),  # SRS_info(organism_droplet_dict[wildcards.organism][0], 'tx'),
		counts = lambda wildcards: well_and_droplet_input(wildcards.organism, wildcards.reference, quant_path, SRS_dict, organism_welltech_dict)
	output:
		cell_info = '{organism}_{reference}_cell_info.tsv',
		matrix = 'quant/{organism}/{reference}_full_sparse_matrix.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/build_sparse_matrix.R {wildcards.organism} {output} {input}
		"""

rule merge_across_references:
	input:
		cell_info = lambda wildcards: expand('{{organism}}_{reference}_cell_info.tsv', reference = ORG_ref(wildcards.organism)),
		matrix = lambda wildcards: expand('quant/{{organism}}/{reference}_full_sparse_matrix.Rdata', reference = ORG_ref(wildcards.organism))
	output:
		'quant/{organism}/full_sparse_matrix.Rdata',
		'{organism}_cell_info.tsv'
	shell:
		"""
		module load R/3.6
		# script analyzes, for macaque, which gene is more detected when using human or macaque ref
		# then creates new matrix blending "best" gene from either human or macaque
		Rscript {git_dir}/src/rebuild_macaque_sparse_matrix.R {wildcards.organism} {output} {input}
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
		#cat {input} | head -n 1 > header
		#mv header {output}
		#grep -hv "^value" {input} >> {output}
		module load R/3.6
		Rscript {git_dir}/src/cat_cell_info.R {input}
		"""
		
rule make_seurat_objs:
	input:
		'cell_info.tsv',
		expand('quant/{organism}/full_sparse_matrix.Rdata', \
			 organism = ['Mus_musculus', 'Macaca_Mulatta', 'Homo_sapiens'])
	output:
		seurat = ('seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__preFilter.seuratV3.Rdata')
	shell:
		"""
		module load R/3.6
		export SCIAD_CONFIG={config_abspath}
		Rscript {git_dir}/src/build_seurat_obj_classic.R \
			{output.seurat} {wildcards.partition} {wildcards.covariate} {wildcards.transform} \
			{wildcards.combination} {wildcards.n_features} {input} 
		"""

rule label_known_cells_with_type:
	input:
		'cell_info.tsv',
		config['srr_sample_file']		
	output:
		'cell_info_labelled.Rdata'
	shell:
		"""
		module load R/3.6
		export SCIAD_CONFIG={config_abspath}
		Rscript {git_dir}/src/label_known_cells.R 
		"""
rule integrate_00:
	input:
		cell_label_info = 'cell_info_labelled.Rdata',
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__preFilter.seuratV3.Rdata',
	output:
		#temp('seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata')
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'
	threads: 2 	
	shell:
		'''
		module load R/3.6
		export SCIAD_CONFIG={config_abspath}
		cmd="Rscript {git_dir}/src/merge_methods.R \
				  {wildcards.method} \
				  {wildcards.transform} \
				  {wildcards.covariate} \
				  {wildcards.dims} \
				  {input.obj} \
				  {output}"
		
		echo $cmd 
		eval $cmd

		'''
	
	# run:
	# 	job = f'module load R/3.6; \
	# 			*Rscript {git_dir}/src/merge_methods.R \
	# 			  {wildcards.method} \
	# 			  {wildcards.transform} \
	# 			  {wildcards.covariate} \
	# 			  {wildcards.dims} \
	# 			  {input} \
	# 			  {output}'
	# 	sp.run("echo \"" +  job + "\"\n", shell = True)
	# 	sp.run(job, shell = True)



# as the wellONLY data has no labelled cells, we need to use the droplet data to transfer labels
def predict_missing_cell_types_input(partition):
	if partition == 'onlyWELL':
		out = ['seurat_obj/{combination}__n_features{n_features}__{transform}__onlyWELL__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata',
				'seurat_obj/{combination}__n_features{n_features}__{transform}__onlyDROPLET__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata']
	else:
		out = ['seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata']
	return(out)

rule predict_missing_cell_types_00:
	input:
		lambda wildcards: predict_missing_cell_types_input(wildcards.partition),
		'cell_info_labelled.Rdata'
	output:
		'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter_cell_info_predictions.Rdata'
	run:
		if wildcards.partition == 'onlyWELL':
			command = 'module load R/3.6; ' \
						'?Rscript {git_dir}/src/transfer_labels__sep_objs.R \
							{input} {transform} {output}'.format(input = ' '.join(input),
																			transform = wildcards.transform,
																			output = output) 
		else:
			command = 'module load R/3.6; ' \
						'Rscript {git_dir}/src/transfer_labels.R \
							{input} {transform} {output}'.format(input = ' '.join(input),
																			transform = wildcards.transform,
																			output = output)
		sp.run("echo \"" +  command + "\"\n", shell = True)
		sp.run(command, shell=True)


rule celltype_predict_VS_xgboost:
	input:
		('seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'),
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umapFilter.Rdata'
	output:
		'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.Rdata',
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umapFilter.predictions.Rdata'		
	shell:
		"""
		module load R/3.6
		export SCIAD_CONDA_DIR={conda_dir}
		export SCIAD_GIT_DIR={git_dir}
		export SCIAD_WORKING_DIR={working_dir}
		Rscript  {git_dir}/src/celltype_predict_wrapper.R {input} {output}
		"""

rule doublet_ID:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'
	output:
		'doublet_calls/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.doublets.Rdata'
	shell:
		"""
		module load R/3.6
		export SCIAD_CONDA_DIR={conda_dir}
		Rscript {git_dir}/src/doublet_ID.R {input} {output}
		"""

rule calculate_umap:
	input:
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'
	output:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	threads: 4
	shell:
		"""
		module load R/3.6
		export SCIAD_CONDA_DIR={conda_dir}
		export SCIAD_GIT_DIR={git_dir}
		Rscript {git_dir}/src/calculate_umap_and_cluster.R \
			{wildcards.method} {wildcards.dims} {wildcards.dist} {wildcards.neighbors} 1 FALSE TRUE {input} {output}
		"""

rule calculate_tsne:
	input:
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'
	output:
		temp('seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__perplexity{perplexity}.tsne.Rdata')
	threads: 4
	shell:
		"""
		module load R/3.6
		export SCIAD_FITSNE_DIR={fi_tsne_dir}
		Rscript {git_dir}/src/calculate_TSNE.R \
			{wildcards.method} {wildcards.dims} {wildcards.perplexity} {input} {output}
		"""

rule calculate_phate:
	input:
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'
	output:
		'phate/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.phate.Rdata'
	threads: 16
	shell:
		"""
		module load R/3.6
		export SCIAD_CONDA_DIR={conda_dir}
		Rscript {git_dir}/src/run_phate.R {input} {output}
		"""
		
rule calculate_cluster:
	input:
		obj = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata'
	output:
		temp('seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.seuratV3.Rdata')
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/calculate_umap_and_cluster.R \
			{wildcards.method} {wildcards.dims} 1 1 {wildcards.knn} TRUE FALSE {input} {output}
		"""

rule extract_umap:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata',
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn7.cluster.Rdata',
		'cell_info_labelled.Rdata'
		#'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter_cell_info_predictions.Rdata'
	output:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/extract_umap.R \
			{input} {output} {wildcards.method}	UMAP
		"""

rule extract_tsne:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__perplexity{perplexity}.tsne.Rdata',
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn7.cluster.seuratV3.Rdata',
		'cell_info_labelled.Rdata'
		#'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter_cell_info_predictions.Rdata'
	output:
		'tsne/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__perplexity{perplexity}.tsne.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/extract_umap.R \
			{input} {output} {wildcards.method} TSNE
		"""

rule extract_cluster:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.seuratV3.Rdata'
	output:
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.graph.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/extract_cluster.R \
			{input} {output}
		"""

rule plot_integration:
	input:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	output:
		'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plot.png'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/big_plots.R UMAP {input} {output}
		"""

rule plot_integration_PREDICT:
	input:
		umap = 'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		celltype_predict = 'umap/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features10000__counts__universe__batch__scVI__dims30__preFilter__mindist0.1__nneighbors100.umapFilter.predictions.Rdata'	
	output:
		'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.big_plotPredictions.png'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/big_plots.R UMAP {input.umap} {output} {input.celltype_predict}
		"""

rule plot_integration_tsne:
	input:
		'tsne/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__perplexity{perplexity}.tsne.Rdata'
	output:
		'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__perplexity{perplexity}.big_tsne_plot.png'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/big_plots.R TSNE {input} {output}
		"""

rule monocle_diff_testing:
	input:
		'monocle_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{knn}.monocle.Rdata'
	output:
		temp('diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__{dims}__{dist}__{knn}__{neighbors}.{piece}.{model}.diff.Rdata')
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/diff_testing_monocle.R {input} 200 {wildcards.piece} {wildcards.model} {output}
		"""

rule monocle_diff_testing_subcluster:
	input:
		'monocle_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{knn}.monocle.Rdata'
	output:
		temp('diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__{dims}__{dist}__{knn}__{neighbors}.G.{cluster_chunk}.diff.Rdata')
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/diff_testing_monocle_subcluster.R {input} 200 {wildcards.cluster_chunk} G {output}
		"""

rule monocle_diff_merge:
	input:
		expand('diff_testing/{{combination}}__{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__{{dims}}__{{dist}}__{{knn}}__{{neighbors}}.{piece}.{{model}}.diff.Rdata', piece = list(range(1,201)))
		#expand('diff_testing/{{combination}}__n_features{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__dims{{dims}}__mindist{{dist}}__knn{{knn}}__nneighbors{{neighbors}}.{piece}pieces.model{{model}}.monocle_diff.Rdata',
		#		piece = list(range(1,201)))
	output:
		'diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__{dims}__{dist}__{knn}__{neighbors}.{model}.diff.coef_table.Rdata',
		#'diff_testing/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__knn{knn}__nneighbors{neighbors}.model{model}.monocle_diff.gene_fits.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/merge_diff_monocle.R {output} {wildcards.model} {input}
		"""

rule monocle_diff_merge_subcluster:
	input:
		expand('diff_testing/{{combination}}__{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__{{dims}}__{{dist}}__{{knn}}__{{neighbors}}.G.{cluster_chunk}.diff.Rdata', cluster_chunk = list(range(1,41)))
	output:
		'diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__{dims}__{dist}__{knn}__{neighbors}.G.SC.diff.coef_table.Rdata',
	shell:
		"""
		module load R/3.6
		?Rscript {git_dir}/src/merge_diff_monocle_subcluster.R {output} G {input}
		"""

rule monocle_marker_test:
	input:
		'monocle_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{knn}.monocle.Rdata'
	output:
		'diff_test/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{knn}__{group}.monocleTopMarker.Rdata'
	threads: 16
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/diff_testing_topMarker_monocle.R {input} {threads} {wildcards.group} {output}
		"""

rule build_monocle_obj:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.Rdata'
	output:
		#'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.monocle_trajectory_labels.png',
		#'plots/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.monocle_trajectory.png',
		'monocle_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{knn}.monocle.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/monocle.R {input} {output}
		"""

rule diff_test_wilcox:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.Rdata',
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata'
	output:
		'diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__{knn}__{neighbors}__{dist}__{group}.sceWilcox.Rdata',
		'diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__{knn}__{neighbors}__{dist}__{group}.sceWilcox_summary.Rdata'
	threads: 4 
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/diff_testing_sce_wilcox.R {input} {wildcards.group} {threads} {output}
		"""

rule perf_metrics:
	input:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist0.3__nneighbors30.umap.Rdata',
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist0.3__nneighbors30.umap.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.Rdata'
	output:
		'perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/perf_metrics.R {input} {output}
		"""

rule make_h5ad_object:
	input:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.Rdata',
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata',
		'cell_info_labelled.Rdata',
		#'predictions/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter_cell_info_predictions.Rdata'
		#'monocle_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{knn}.monocle.Rdata'
	output:
		temp('anndata/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}.h5ad')
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/seurat_to_h5ad.R {input} {output}
		"""

rule scIB_stats:
	input:
		 'anndata/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}.h5ad'
	output:
		'scIB_stats/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}___stats.csv'
	shell:
		"""
		module load R/3.6
		export SCIAD_CONDA_DIR={conda_dir}
		export SCIAD_GIT_DIR={git_dir}
		Rscript {git_dir}/src/scIB_stats.R {wildcards.method} {input} {wildcards.dims} {output}
		"""

rule doublet_filtering:
	input:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		'doublet_calls/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.doublets.Rdata'
	output:
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umapFilter.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/doublet_filtering.R {input} {output}
		"""

rule make_sqlite:
	input:
		seurat = 'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umap.Rdata',
		meta = 'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umapFilter.predictions.Rdata',
		cluster = 'cluster/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.cluster.Rdata',
		well_data =  'well_data_seurat_obj_labelled.Rdata',
		phate = 'phate/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.phate.Rdata'
	params:
		'site/anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite'
	output:
		'site/anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite'
	shell:
		"""
		module load R/3.6
		export SCIAD_GIT_DIR={git_dir}
		Rscript {git_dir}/src/make_sqlite.R \
			{input.seurat} \
			{input.meta} \
			{input.cluster} \
			{input.well_data} \
			{input.phate} \
			{output} \
			{wildcards.correction}
		#pigz -p 32 {params}
		"""

rule pseudoBulk_DGE_buildObj:
	input:
		'seurat_obj/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter.seuratV3.Rdata',
		'umap/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}.umapFilter.predictions.Rdata',
		'well_data_seurat_obj_labelled.Rdata'
	output:
		'pseudoBulk_DGE/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{pseudoTest}__edgeR_obj.Rdata'
	threads: 6
	shell:
		"""
		module load R/4.0
		export SCIAD_GIT_DIR={git_dir}
		Rscript {git_dir}/src/pseudoBulk_buildObj.R {input} {wildcards.pseudoTest} {output}		
		"""

rule pseudoBulk_DGE_difftest:
	input:
		'pseudoBulk_DGE/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{pseudoTest}__edgeR_obj.Rdata'
	output:
		'pseudoBulk_DGE/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__{pseudoTest}__{piece}.Rdata'
	shell:
		"""
		module load R/4.0
		export SCIAD_GIT_DIR={git_dir}
		Rscript {git_dir}/src/pseudoBulk_diff_testing.R {input} {wildcards.pseudoTest} {wildcards.piece} {output}		
		"""
	
rule merge_pseudoBulk_ABC:
	input:
		pseudoBulk = expand('pseudoBulk_DGE/{{combination}}__n_features{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__dims{{dims}}__preFilter__mindist{{dist}}__nneighbors{{neighbors}}__{pseudoTest}__{piece}.Rdata', pseudoTest = ['A1','A2','A3','B1','B2','B3','C1','C3'], piece = range(1,26)),
		        pseudoBulk_well = expand('pseudoBulk_DGE/{{combination}}__n_features{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__dims{{dims}}__preFilter__mindist{{dist}}__nneighbors{{neighbors}}__{pseudoTest}__{piece}.Rdata', pseudoTest = ['Cw1','Cw2','Cw3'], piece = range(1,2))
	output:
		'pseudoBulk_DGE/merge/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__ABC.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/pseudoBulk_merge.R {input} {output}
		"""

rule merge_pseudoBulk_C2:
	input:
		pseudoBulk = expand('pseudoBulk_DGE/{{combination}}__n_features{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__dims{{dims}}__preFilter__mindist{{dist}}__nneighbors{{neighbors}}__{pseudoTest}__{piece}.Rdata', pseudoTest = ['C2'], piece = range(1,501))
	output:
		'pseudoBulk_DGE/merge/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__C2.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/pseudoBulk_merge.R {input} {output}
		"""


rule sqlite_add_tables:
	input:
		sqlite = 'site/anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite',
		#diff_wilcox = expand('diff_testing/{{combination}}__{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__dims{{dims}}__{{knn}}__{{neighbors}}__{{dist}}__{group}.sceWilcox.Rdata', \
		#			group = ['subcluster', 'cluster','CellType_predict','CellType']),
		#diff_glm = expand('diff_testing/{{combination}}__{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__{{dims}}__{{dist}}__{{knn}}__{{neighbors}}.{model}.diff.coef_table.Rdata', \
		#			model = ['A','C', 'E']),
		#diff_glm_subcluster = 'diff_testing/{combination}__{n_features}__{transform}__{partition}__{covariate}__{method}__{dims}__{dist}__{knn}__{neighbors}.G.SC.diff.coef_table.Rdata',
		#marker_monocle = expand('diff_test/{{combination}}__n_features{{n_features}}__{{transform}}__{{partition}}__{{covariate}}__{{method}}__dims{{dims}}__preFilter__mindist{{dist}}__nneighbors{{neighbors}}__{{knn}}__{group}.monocleTopMarker.Rdata', \
		#			group = ['seuratCluster','CellType_predict']),
		ABC = 'pseudoBulk_DGE/merge/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__ABC.Rdata',
		C2 = 'pseudoBulk_DGE/merge/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__C2.Rdata',
		doublet = 'doublet_calls/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}.doublets.Rdata'
	output:
		uncompressed = 'site/MOARTABLES__anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite',
		compressed = 'site/MOARTABLES__anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite.gz'
	params:
		inp = 'site/anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite',
		uncompressed = 'site/MOARTABLES__anthology_limma{correction}___{combination}-{n_features}-{transform}-{partition}-{covariate}-{method}-{dims}-{dist}-{neighbors}-{knn}.sqlite'
	threads: 16
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/sqlite_add_diff_tables.R {params.inp} \
			{input.ABC} \
			{input.C2} \
			{input.doublet}
		pigz -c -p {threads} {params.inp} > {output.compressed}
		mv {input.sqlite} {params.uncompressed}
		"""

rule make_PLAE_objs:
	input:
		'site/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite',
		'seurat_obj/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features1000__counts__raw__batch__preFilter.seuratV3.Rdata'
	output:
		'site/scEiaD_droplet_seurat_v3.Rdata',
		'site/scEiaD_well_seurat_v3.Rdata',
		'site/counts_unfiltered.Rdata',
		'site/counts.Rdata',
		'site/cpm.Rdata'
	shell:
		"""
		module load R/3.6
		Rscript {git_dir}/src/output_objs_for_plae.R

		Rscript {git_dir}/src/seurat_to_h5ad_core.R site/scEiaD_droplet_seurat_v3.Rdata scEiaD_droplet site/scEiaD_droplet_anndata.h5ad
		Rscript {git_dir}/src/seurat_to_h5ad_core.R site/scEiaD_well_seurat_v3.Rdata scEiaD_well site/scEiaD_well_anndata.h5ad

		Rscript {git_dir}/src/pseudoBulk_output_for_site.R pseudoBulk_DGE/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15__A1__edgeR_obj.Rdata site/pseudoBulk_celltype.tsv.gz
		Rscript {git_dir}/src/pseudoBulk_output_for_site.R pseudoBulk_DGE/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15__B1__edgeR_obj.Rdata site/pseudoBulk_celltypePredict.tsv.gz
		Rscript {git_dir}/src/pseudoBulk_output_for_site.R pseudoBulk_DGE/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15__C1__edgeR_obj.Rdata site/pseudoBulk_clusterDroplet.tsv.gz
		Rscript {git_dir}/src/pseudoBulk_output_for_site.R pseudoBulk_DGE/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_features5000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.1__nneighbors15__Cw1__edgeR_obj.Rdata site/pseudoBulk_clusterWell.tsv.gz

		Rscript {git_dir}/src/build_QC_stats.R
		"""


if config['subset_clustering'] == 'False':	
	rule merge_stats:
		input:	
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['counts'], \
					method = ['scVI'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['universe'], \
					n_features = [2000, 5000, 10000], \
					covariate = ['batch'], \
					dims = [30,50,100],
					knn= [7]),
			expand('scIB_stats/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}___stats.csv', \	
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					n_features = [2000, 5000, 10000], \
					transform = ['counts'], \
					partition = [ 'universe'], \
					covariate = ['batch'], \
					method = ['scVI'], \
					dims = [30,50,100], \
					dist = [0.1], \
					neighbors = [100], \
					knn = [7]),
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['counts'], \
					method = ['scVI'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['TabulaDroplet'], \
					n_features = [1000, 2000, 5000, 10000], \
					covariate = ['batch'], \
					dims = [4,6,8,10,20,30,50,100],
					knn = [0.2,0.4,0.6,5, 7, 10]),
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['libSize','sqrt','scran', 'standard'], \
					method = ['bbknn','insct',  'magic', 'scanorama', 'harmony', 'fastMNN', 'combat',  'none'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['TabulaDroplet'], \
					n_features = [2000], \
					covariate = ['batch'], \
					dims = [8, 30],
					knn = [7]),
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['libSize','sqrt','scran', 'standard'], \
					method = ['bbknn','insct','magic', 'scanorama', 'harmony', 'fastMNN', 'combat', 'liger', 'none'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['onlyWELL'], \
					n_features = [500, 1000, 2000], \
					covariate = ['batch'], \
					dims = [4, 8, 30],
					knn = [7]),
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['counts'], \
					method = ['scVI'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['onlyWELL'], \
					n_features = [500, 1000, 2000], \
					covariate = ['batch'], \
					dims = [4, 8, 30],
					knn = [7]),
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['SCT'], \
					method = ['CCA'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['onlyWELL'], \
					n_features = [500, 1000, 2000], \
					covariate = ['batch'], \
					dims = [4, 8, 30],
					knn = [7]),
			expand('scIB_stats/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}___stats.csv', \	
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					n_features = [1000, 2000, 5000, 10000], \
					transform = ['counts'], \
					partition = [ 'TabulaDroplet'], \
					covariate = ['batch'], \
					method = ['scVI'], \
					dims = [4,6,8,10,20,30,50,100], \
					dist = [0.1], \
					neighbors = [500], \
					knn = [0.2,0.4,0.6,5,7,10]),
			expand('scIB_stats/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}___stats.csv', \
					transform = ['libSize','sqrt','scran', 'standard'], \
					method = ['bbknn','insct', 'magic', 'scanorama', 'harmony', 'fastMNN', 'combat', 'CCA', 'none'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['TabulaDroplet'], \
					n_features = [2000], \
					covariate = ['batch'], \
					dims = [8, 30], \
					knn = [7], \
					dist = [0.3], \
					neighbors = [30])
		output:
			'merged_stats_2020_07_06.Rdata'
		shell:
			"""
			module load R/3.6
			Rscript {git_dir}/src/optimal_params.R
			"""

else: 
	rule merge_stats:
		input:	
			expand('perf_metrics/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__knn{knn}.Rdata', \
					transform = ['counts'], \
					method = ['scVI'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['cones', 'hc', 'rgc', 'amacrine', 'mullerglia', 'bipolar', 'rods'], \
					n_features = [ 1000, 2000], \
					covariate = ['batch'], \
					dims = [4,6,10,30],
					knn = [7, 10, 15]),
			expand('scIB_stats/{combination}__n_features{n_features}__{transform}__{partition}__{covariate}__{method}__dims{dims}__preFilter__mindist{dist}__nneighbors{neighbors}__knn{knn}___stats.csv', \
					transform = ['counts'], \
					method = ['scVI'], \
					combination = ['Mus_musculus_Macaca_fascicularis_Homo_sapiens'], \
					partition = ['cones', 'hc', 'rgc', 'amacrine', 'mullerglia', 'bipolar', 'rods'], \
					n_features = [ 1000, 2000], \
					covariate = ['batch'], \
					dims = [4,6,10,30], \
					knn = [7, 10, 15], \
					dist = [0.1], \
					neighbors = [30])
		output:
			'merged_stats_subsetClustering.Rdata'
		shell:
			"""
			module load R/3.6
			Rscript {git_dir}/src/optimal_params.R
			"""
