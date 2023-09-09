import scanpy as sc

big_h5ad = 'site/scEiaD_all_anndata.h5ad'

adata = sc.read_h5ad(big_h5ad)

well_counts = adata[adata.obs['TechType'].isin(['SMARTSeq_v2'])].obs[['batch']].value_counts()
well_batches = list(well_counts[well_counts > 150].reset_index()['batch'])

with open(r'well_batches.txt', 'w') as fp:
	for item in well_batches:
		fp.write("%s\n" % item)


droplet_counts = adata[~adata.obs['TechType'].isin(['SMARTSeq_v2'])].obs[['sample_accession']].value_counts()
droplet_samples = list(droplet_counts[droplet_counts > 150].reset_index()['sample_accession'])


with open(r'droplet_counts.txt', 'w') as dc:
	for item in droplet_samples:
		if item == 'NA':
			continue
		dc.write("%s\n" % item)
