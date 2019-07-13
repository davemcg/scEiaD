import glob as glob
import subprocess as sp

srr_sample_file = config['srr_sample_file']

# builds dictionary of dictionaries where first dict key is SRS 
# and second dict key are SRS properties
SRS_dict={}
with open(srr_sample_file) as file:
	for line in file:
		info = line.strip('\n').split('\t')
		if info[0] == 'sample_accession':
			continue
		SRS = info[0]
		if SRS not in SRS_dict:
			SRS_dict[SRS]={'SRR': [info[1]],
					  'paired':True if info[2]=='PAIRED' else False, 
				      'organism':info[3],
		              'tech':info[4],
					  'UMI':True if info[5]=='YES' else False}
		else:
			runs = SRS_dict[SRS]['SRR']
			runs.append(info[1])
			SRS_dict[SRS]['SRR'] = runs

def lookup_run_from_SRS(SRS):
	#i = '0'
	#print(SRS)
	#if '_' in SRS:
	#	i= '1' if SRS[-1]=='1' else '2'# check L/R file
	#	SRS=SRS[:-2]
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
	print(out)
	return(out)

def SRS_kallisto_idx(SRS):
	organism = SRS_dict[SRS]['organism']
	if organism.lower() == 'mus musculus':
		idx = 'references/kallisto_idx/gencode.vM22.pc_translations.fa.gz.idx'
	elif organism.lower() == 'homo sapiens':
		idx = 'references/kallisto_idx/gencode.v31.pc_transcripts.fa.gz.idx'
	elif organism.lower() == 'macaca fascicularis':
		idx = 'references/kallisto_idx/GCF_000364345.1_Macaca_fascicularis_5.0_protein.faa.gz.idx'
	else:
		print(SRS + ' ' + organism + " NO SPECIES MATCH!")
	return(idx)

SRS_UMI_samples = []
for SRS in SRS_dict.keys():
	if SRS_dict[SRS]['UMI']:
		SRS_UMI_samples.append(SRS)

wildcard_constraints:
	SRS = '|'.join(SRS_UMI_samples)

rule all:
	input:
		expand('fastq/{SRS}.fastq.gz', SRS = SRS_UMI_samples)
		#expand('quant/{SRS}/output.bus', SRS = SRS_UMI_samples)
		#expand('fastq/{SRS}.fastq.gz', SRS = SRS_UMI)
		#expand(SRR_FASTQ)
		#expand('{srr}{layout}.fastq.gz', srr = SRR_IDs, layout = library_layout(wildcards.srr))

# aggregate multiple SRR runs into one SRS sample
rule SRR_to_SRS:
	input:
		lambda wildcards: lookup_run_from_SRS(wildcards.SRS)
	output:
		'fastq/{SRS}.fastq.gz'
	run:	
		SRR_files = lookup_run_from_SRS(SRS)
		i = '1' if '_' in SRS and SRS[-1]=='1' else '2'# which strand
		SRS = SRS[:-2] if '_' in SRS else SRS
		for SRR in SRR_files:
			if SRS_dict[SRS]['paired']:
				sp.run('cat {SRR} >> fastq_files/{SRS}_{i}.fastq.gz'.format(SRR=SRR, i=i, SRS=SRS), shell=True)
			else:
				sp.run('cat {SRR} >> fastq_files/{SRS}.fastq.gz'.format(SRR=SRR, SRS=SRS), shell=True)

# mouse, human, macaque fasta and gtf
rule download_references:
	output:
		mouse_fasta = 'references/gencode.vM22.pc_translations.fa.gz',
		macaque_fasta = 'references/GCF_000364345.1_Macaca_fascicularis_5.0_protein.faa.gz',
		human_fasta = 'references/gencode.v31.pc_transcripts.fa.gz'
	shell:
		"""
		mkdir -p references
		wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.pc_translations.fa.gz  
		wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_protein.faa.gz'
		wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.pc_transcripts.fa.gz
		mv *fa*gz references/
		"""

# need to make mouse, human, macaque
rule kallisto_index:
	input:
		'references/{fasta}'
	output:
		'references/kallisto_idx/{fasta}.idx'
	conda:
		'envs/kallisto.yaml'
	shell:
		"""
		kallisto index {input} -i {output}
		"""

# use the gtf to make the bustool count tx file
# tsv three column
# transcript_ID gene_ID gene_name
# ENST ENSG gene_name
#rule tx_gene_mapping:

# this does the pseudoalignment
rule kallisto_bus:
	input:
		fastq = ancient(lambda wildcards: ['fastq/{}_1.fastq.gz'.format(wildcards.SRS), \
											'fastq/{}_2.fastq.gz'.format(wildcards.SRS)] \
											if SRS_dict[wildcards.SRS]['paired'] \
											else 'fastq/{}.fastq.gz'.format(wildcards.SRS)),
		idx = lambda wildcards: SRS_kallisto_idx(wildcards.SRS)
	conda:
		'envs/kallisto.yaml'
	output:
		'quant/{SRS}/output.bus'
	params:
		tech = lambda wildcards: SRS_dict[wildcards.SRS]['tech']
	shell:
		"""
		kallisto bus {input.fastq} -x {params.tech} -i {input.idx} -o quant/{wildcards.SRS}
		"""
	
# sorting required for whitelist creation and correction
# make these temp files
#rule bustools_sort:
#	input:
#	output:
#	shell:
#		"""
#		/home/mcgaugheyd/git/bustools/build/src/./bustools
#		"""

# find barcodes, correct barcodes
# make these temp files
#rule bustools_whitelist_correct:

# does the gene level quant
#rule bustools_count:

#rule profit:
