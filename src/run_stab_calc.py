#%%
import pandas as pd
import numpy as np 
import glob
import os 
import re
import pickle
from clustereval import stability
from multiprocessing import Pool
os.chdir('/data/swamyvs/scEiad_subcelltype')
#%%
NCORES = 64

exp_files = glob.glob("clustering_out/*/*.csv.gz", recursive=True)


exp_meta_df = pd.DataFrame([re.split("_|-|/", x.split('out/')[1])[::2] for x in exp_files],
                           columns=['clustering_type', 'experiment', 'dist', 'knn']).assign(path=exp_files)
#%%
exp_dfs = [pd.read_csv(i).assign(
    labels=lambda x: 'clu_' + x.labels.astype(str)) for i in exp_files]

print("Calculating H_tot\n")
exp_meta_df = exp_meta_df.assign(H_tot=[stability.entropy(i) for i in exp_dfs])



exp_stability_gen = (
    (exp_dfs[i], 
     exp_meta_df.drop(axis=1, index=i), exp_dfs[0:i] +
     exp_dfs[(i+1):], 
     '-'.join(exp_meta_df.iloc[i, :4].to_list()))
    for i in range(len(exp_dfs))
)

print("Running Stablity calc")
pool = Pool(NCORES)
exp_stability_dfs = pool.map(stability.calc_stability, exp_stability_gen)
exp_stability_dfs = [i for i in exp_stability_dfs]



# %%
with open("testing/exp_stab_df_list_full_try4_noperturb.pck", 'wb+') as ofl:
    pickle.dump(exp_stability_dfs, ofl)

print(exp_stability_dfs)

# %%
