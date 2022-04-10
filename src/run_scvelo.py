import scanpy as sc
import scvelo as scv
import argparse
parser = argparse.ArgumentParser(description='Run scVelo and save new adata object')
parser.add_argument('adata_file', help='Input adata file')
parser.add_argument('adata_out', help='Output adata prefix (will output PREFIX.SAMPLE.h5ad)')
parser.add_argument('--redo_pca', help='Use exising PCA, (Y)es or (N)o')
parser.add_argument('--move_scvi_to_pca', help='Copy scVI latent dims into pca slot, (Y)es or (N)o')

args = vars(parser.parse_args())
print(len(args))
print(args)

redo_pca = False
move_scvi_to_pca = False

adata_file = args['adata_file']
adata_out = args['adata_out']
if args['redo_pca'] == 'Y':
	print("Redo PCA")
	redo_pca = True
if args['move_scvi_to_pca']  == 'Y':
	print("Use scVI as PCA")
	move_scvi_to_pca = True

def run_scVelo(adata, adata_out, redo_pca, move_scvi_to_pca, sample_name):
	print(adata)
	n_pcs = 20
	if redo_pca:
		print("\nRunning HVG and PCA\n")
		scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=10000)
		sc.tl.pca(adata)
	if move_scvi_to_pca:
		adata.obsm['X_pca'] = adata.obsm['X_scVI']
		n_pcs = adata.obsm['X_pca'].shape[1]
	sc.pp.neighbors(adata, n_pcs=n_pcs)
	# sc.tl.umap(adata)
	scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
	print("Recover dynamics")
	scv.tl.recover_dynamics(adata, n_jobs=12)
	scv.tl.velocity(adata, mode="dynamical")
	scv.tl.velocity_graph(adata, n_jobs=4)
	# scv.tl.velocity_embedding(adata, basis="umap")
	# https://github.com/theislab/scvelo/issues/255#issuecomment-739995301
	adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
	adata.write_h5ad(adata_out + '.' + sample + '.h5ad')
	return(adata)

def scVelo_plotting(adata_path, prefix):
	adata = sc.read_h5ad(adata_path)
	#scv.tl.velocity_embedding(adata, basis="umap")
	scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=4,  color='CellType_predict', save= prefix + '.embedding_stream.png', show = False, dpi = 300)
	#scv.pl.velocity(adata, ['SALL1'], show = False, save = prefix + '.SALL1.png', dpi = 300, color ='seurat_clusters',  palette=color_palette)
	scv.tl.rank_velocity_genes(adata, groupby='CellType_predict', min_corr=.3)
	df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
	df.to_csv(prefix + '.rank_velocity.csv')
	adata.var.to_csv(prefix + '.var.csv')
	scv.tl.paga(adata, groups='CellType_predict')
	scv.pl.paga(adata, basis='umap', size=4, alpha=.3, color = 'CellType_predict',
		min_edge_width=2, node_size_scale=1.5, show = False,
		dpi = 300, save = prefix + '.PAGA.png')

adata = sc.read_h5ad(adata_file)
samples = [x.split('_') for x in adata.obs['Barcode'].to_numpy()]
samples = [x[1] for x in samples if len(x) == 2]
uniq_samples = list(set(samples))

for sample in uniq_samples:
	adata_sample = adata[[sample in x for x in adata.obs['Barcode']], :].copy()
	run_scVelo(adata_sample, adata_out, redo_pca, move_scvi_to_pca, sample)

for sample in uniq_samples:
	scVelo_plotting('scvelo/' + adata_out + '.' + sample + '.h5ad', 'scvelo/' + sample)
