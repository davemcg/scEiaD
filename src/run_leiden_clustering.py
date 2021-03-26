

import numpy as np
import pandas as pd
import time
import os 
import sys
import leidenalg
from joblib import Parallel, delayed
from clustereval import cluster 


def leiden_clustering(reduction, res, k, i):
    clu_obj = cluster.Cluster(data=reduction, knn=k,  nthreads=1)
    clu_obj.buildNeighborGraph(nn_space='l2', ef_construction=150,
                            local_pruning=True, global_pruning=True, jac_std_global='median')
    #clu_obj.run_perturbation()
    labels = clu_obj.run_leiden(
        vertex_partition_method=leidenalg.RBConfigurationVertexPartition,
        n_iter=5,
        resolution=res,
        jac_weighted_edges='weight'
    )
    labels_corrected = clu_obj.merge_singletons(labels, 10)
    outdf = pd.DataFrame(
        {"Barcode": list(reduction.index), 'labels': labels_corrected})
    outdf.to_csv(
        f'clustering_out/leiden/exp-{str(i)}_resolution-{str(res)}_knn-{str(k)}_.csv.gz', index=False)
    return


def run_exp(reduction_file, dist):
    reduction = pd.read_csv(reduction_file, index_col=0)
    knn_uniform = list(range(15, 101))
    Parallel(n_jobs=-1)(
        delayed(leiden_clustering)(reduction, dist, k, i)
        for i, k in enumerate(knn_uniform))


# %%


run_exp('testing/amacrin_mm_scvi_dim.csv', float(sys.argv[1]))

# %%
