#%%
import pandas as pd
import numpy as np 
import glob
import os 
import re
from joblib import Parallel, delayed
import pickle
from concurrent.futures import ProcessPoolExecutor
os.chdir('/data/swamyvs/scEiad_subcelltype')

#%%


exp_files = glob.glob("parc_exp/*.csv.gz")


exp_meta_df = pd.DataFrame([ re.split("_|-", x)[2:7: 2] for x in exp_files], columns = ['experiment', 'dist', 'knn']).assign(path = exp_files)




def entropy(exp): # both are dfs with two columsn, Barcode,cluster
    # calc H_tot
    entropy = (exp
             .groupby("labels")
             .count()
             .reset_index(drop=True)
             .assign(prop=lambda x: x/exp.shape[0],
                     H=lambda x: x['prop'] * np.log(x['prop']))
             ['H'].sum()*-1)
    return entropy
    

def H_k(ref, exp):
    exp = exp[exp.Barcode.isin(ref['Barcode']) ]
    if exp.shape[0] == 0:
        return 0
    else:
        h_k = entropy(exp)
        return h_k

# %%

exp_dfs = [pd.read_csv(i) for i in exp_files]

# %%
exp_meta_df = exp_meta_df.assign(H_tot = [entropy(i) for i in exp_dfs ])

# %%

sanes_lab_df = pd.read_csv('testing/sanes_bc_lab.csv')

# %%
# this is the sanes level scores 


def calc_stability(tup):
    ref = tup[0]
    meta_df = tup[1]
    exp_df_list = tup[2]
    try:
        H_k_scores = np.asarray([
                                    [H_k(group[1], exp) for exp in exp_df_list] / meta_df['H_tot']
                                    for group in ref.groupby("labels") 
                                ] ).sum(axis = 1) 
        H_k_scores = 1 - (H_k_scores / len(exp_df_list)) 
        clusters = [group[0] for group in ref.groupby("labels") ]
        return pd.DataFrame.from_dict({'clusters' : clusters , 'H_k_scores' : H_k_scores})
    except:
        return None

#sanes_clustering_stability = calc_stability(sanes_lab_df, exp_meta_df, exp_dfs)


#%%
NCORES = 16



exp_stability_gen = (
    (exp_dfs[i],exp_meta_df.drop(axis = 1, index=i), exp_dfs[0:i] + exp_dfs[(i+1):])
    for i in range(len(exp_dfs))
)

with ProcessPoolExecutor(max_workers=NCORES) as proc_exec:
    exp_stability_dfs = proc_exec.map(calc_stability, exp_stability_gen)
    exp_stability_dfs = [i for i in exp_stability_dfs]



# %%
with open("testing/exp_stab_df_list.pck", 'wb+') as ofl:
    pickle.dump(exp_stability_dfs, ofl)

# %%


