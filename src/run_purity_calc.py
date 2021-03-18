#%%
import pandas as pd
import numpy as np 
import glob
import os 
import re
import pickle
from clustereval import purity
from multiprocessing import Pool
os.chdir('/data/swamyvs/scEiad_subcelltype')

#%%

exp_files = glob.glob("clustering_out/*/*.csv.gz", recursive=True)
exp_meta_df = pd.DataFrame([re.split("_|-|/", x.split('out/')[1])[::2] for x in exp_files],
                           columns=['clustering_type', 'experiment', 'dist', 'knn']).assign(path=exp_files)
exp_dfs = [pd.read_csv(i).assign(
    labels=lambda x: 'clu_' + x.labels.astype(str)) for i in exp_files]

sanes_lab_df = pd.read_csv('testing/sanes_bc_lab.csv')
exp_cluster_sets = [dict(purity.df2DictSets(df)) for df in exp_dfs]
sanes_cluster_sets = purity.df2DictSets(sanes_lab_df)

# %%
exp_purity_gen = (
    (exp_cluster_sets[i],
     exp_cluster_sets[0:i] + exp_cluster_sets[(i+1):],
     '-'.join(exp_meta_df.iloc[i, :4].to_list())
     )
    for i in range(len(exp_cluster_sets))
)

NCORES=48
pool = Pool(NCORES)
res = pool.map(purity.run_purity_calc, exp_purity_gen)
res= [i for i in res]
with open("testing/exp_purity_df_list_full_try3.pck", 'wb+') as ofl:
    pickle.dump(res, ofl)
# %%
