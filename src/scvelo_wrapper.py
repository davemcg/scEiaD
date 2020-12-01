import scvelo as scv
import scanpy as sc
import pandas as pd
def run_velocity(adata,embedding_key, dimred_key):
    adata.layers['spliced'] = adata.X.copy()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    scv.pp.moments(adata, use_rep=embedding_key, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_embedding(adata, basis=dimred_key, autoscale=False)
    vkey=f'velocity_{dimred_key}'
    vdf = pd.DataFrame(adata.obsm[vkey]).assign(pt=list(adata.obs['velocity_pseudotime']))
    return vdf

