#%%
import pandas as pd
import numpy as np 
import glob
import os 
import re
import pickle
from multiprocessing import Pool
os.chdir('/data/swamyvs/scEiad_subcelltype')

#%%

def df2DictSets(df):
    dict_list = df.groupby("labels").apply(
        lambda x: {x['labels'].to_list()[0]: x['Barcode'].to_list()})

    outdict = {key:set(value) for d in dict_list for key,value in d.items() }
    return outdict


exp_files = glob.glob("parc_exp/*.csv.gz")


exp_meta_df = pd.DataFrame([ re.split("_|-", x)[2:7: 2] for x in exp_files], columns = ['experiment', 'dist', 'knn']).assign(path = exp_files)
exp_dfs = [pd.read_csv(i).assign(labels = lambda x: 'clu_' + x.labels.astype(str)) for i in exp_files]

sanes_lab_df = pd.read_csv('testing/sanes_bc_lab.csv')


# %%
exp_cluster_sets = [dict(df2DictSets(df) ) for df in exp_dfs]
sanes_cluster_sets = df2DictSets(sanes_lab_df)
# %%

def purity_k(ref_barcodes, ref_name, query_dict):
    max_overlap = 0
    max_overlap_cluster = "NA"
    for query_cluster in query_dict.keys():
        overlap = ref_barcodes.intersection(query_dict[query_cluster])
        if len(overlap ) > max_overlap:
            max_overlap =len(overlap )
            max_overlap_cluster = query_cluster
    if max_overlap_cluster == "NA":
        print(f"No overlap for {ref_name} in current query")
        return 0
    else:
        purity = max_overlap / len(query_dict[max_overlap_cluster])
        return purity



def run_purity_calc(tup):
    try:
        ref_dict=tup[0]
        query_dict_list=tup[1]
        all_ref_cluster_purities = []
        for ref_cluster in ref_dict.keys():
            ref_cluster_purities = [None] * len(query_dict_list)
            i=0
            for query_dict in query_dict_list:
                ref_cluster_purities[i] = purity_k(ref_dict[ref_cluster], ref_cluster, query_dict)
                i+=1
            all_ref_cluster_purities.append(ref_cluster_purities)
        avg_cluster_purities = pd.DataFrame(
            all_ref_cluster_purities).sum(axis=1) / len(query_dict_list)
        return pd.DataFrame.from_dict({"cluster": list(ref_dict.keys()), "purity": avg_cluster_purities})
    except:
        print("FAIL")
        return None

    

# %%
exp_purity_gen = (
    (exp_cluster_sets[i],
     exp_cluster_sets[0:i] + exp_cluster_sets[(i+1):])
    for i in range(len(exp_cluster_sets))
)

NCORES=48
pool = Pool(NCORES)
res = pool.map(run_purity_calc, exp_purity_gen)
res= [i for i in res]
with open("testing/exp_purity_df_list_full_try2.pck", 'wb+') as ofl:
    pickle.dump(res, ofl)



# %%
res[1]
# %%

# %%
res
# %%
