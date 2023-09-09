import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import argparse
from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser(
                    prog='csv input to seacell output',
                    description='Takes csv and outputs seacell metacells',
                    epilog='David McGaughey 2023 CC0 US Government')

parser.add_argument('h5ad')
parser.add_argument('sample_name')
parser.add_argument('output_obs_name')
parser.add_argument('output_seacell_name')
parser.add_argument('filter_on_col')
 
args = vars(parser.parse_args())
print(args)

h5ad = args['h5ad']
sample = args['sample_name']
obs_file = args['output_obs_name']
seacell_file = args['output_seacell_name']
filter_on_col = args['filter_on_col']

adata = sc.read_h5ad(h5ad)

# filter down to one study or covariate
adata = adata[adata.obs[filter_on_col] == sample]
adata.X = csr_matrix(adata.X)
print(adata)

raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500)

# Compute principal components -
sc.tl.pca(adata, n_comps=20, use_highly_variable=True)
# many errors around this number of cells
# if so, then split ad into pieces
n = 8000 #max num of cells in an adata
num_splits = int(np.ceil(len(adata) / n))
# if num_splits > 1 then tweak n to make equal sized objects
if num_splits > 1:
	n = int(np.ceil(len(adata) / num_splits))

# shuffle to randomize the cells
#adata = adata[np.random.permutation(len(adata))]
# Split the AnnData object into sub-AnnData objects
splits = []
for i in range(num_splits):
    start_idx = i * int(n)
    end_idx = min(start_idx + int(n), len(adata))
    splits.append(adata[start_idx:end_idx])

obs_list = []
seacell_list = []

ad_counter = 0
for ad in splits:
	## Core parameters
	## want ~ 75 cells per meta
	## but having fewer than 10 seacells tends to throw errors
	n_SEACells = max(10, int(ad.shape[0] / 75))
	build_kernel_on = 'X_pca'
	## Additional parameters
	n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

	model = SEACells.core.SEACells(ad,
                    build_kernel_on=build_kernel_on,
                    n_SEACells=n_SEACells,
                    n_waypoint_eigs=n_waypoint_eigs,
                    convergence_epsilon = 1e-5)

	model.construct_kernel_matrix()

	M = model.kernel_matrix
	# Initialize archetypes
	model.initialize_archetypes()
	# fit model
	model.fit(min_iter=5, max_iter=200)
	#model.plot_convergence()
	labels,weights = model.get_soft_assignments()
	meta_aggr = SEACells.core.summarize_by_soft_SEACell(ad, model.A_,summarize_layer='raw', minimum_weight=0.05)
	meta_aggr_df = pd.DataFrame(meta_aggr.X.toarray())
	meta_aggr_df.columns = meta_aggr.var_names
	meta_aggr_df.index = sample + '__' + str(ad_counter) + '__' +  meta_aggr_df.index.astype('str')
	
	ad.obs[['counter']] = ad_counter
	ad.obs[['id']] = sample + '__' + str(ad_counter) + '__' + ad.obs.SEACell.str.extract('(\d+)')
	obs_list.append(pd.DataFrame(ad.obs))
	seacell_list.append(meta_aggr_df)
	
	ad_counter += 1



pd.concat(obs_list).to_csv(obs_file)
pd.concat(seacell_list).to_csv(seacell_file)
