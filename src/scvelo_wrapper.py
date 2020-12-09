import scvelo as scv
import scanpy as sc
import pandas as pd
def run_velocity(adata,embedding_key, dkey, vrank_gb):
    #adata.write_h5ad('testing/tabulaDroplet_scvel0.h5ad')
    adata.layers['spliced'] = adata.X.copy()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    ekey=f'X_{embedding_key}'
    try:
        scv.pp.moments(adata, use_rep=ekey, n_neighbors=30)
    except:
        print('saving neighbor graph')
        scv.pp.neighbors(
            adata=adata,
            n_neighbors=30,
            use_rep=embedding_key,
            use_highly_variable=True,
            n_pcs=None,
            method='umap'
        )
        adata.write_h5ad('testing/tabulaDroplet_scvel0_neighbors.h5ad')
    
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_embedding(adata, basis=dkey, autoscale=False)
    scv.tl.rank_velocity_genes(adata, groupby=vrank_gb, min_corr=.3)
    n_uns = { 'rank_velocity_genes_names' : pd.DataFrame(adata.uns['rank_velocity_genes']['names']), 
          'rank_velocty_genes_scores' : pd.DataFrame(adata.uns['rank_velocity_genes']['scores'])
        }
    adata_for_transport = adata.copy()
    del adata_for_transport.obsp
    adata_for_transport.uns = n_uns
    return adata_for_transport

