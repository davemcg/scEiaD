#!/usr/bin/env python3
#%% 
import os
import sys
import json
import pickle
import re 
from snakemake.utils import read_job_properties
#%%
cluster_json_file =sys.argv[1]
jobscript = sys.argv[2]
custom_config_rules = ['integrate_00', 'calculate_umap', 'calculate_cluster', 'perf_metrics']
#%%

#%%
with open(cluster_json_file) as j:
    cluster_json = json.load(j)
#%%

params = cluster_json['__default__']
job_properties = read_job_properties(jobscript)
params['cpus-per-task'] = job_properties['threads']# we are setting job properties within the snakefile
#%%

# add support for groups so we can bundle multiple rules together and stop dropping our job priorities
if 'rule' in job_properties:
    rule = job_properties['rule']
else:
    rule = job_properties['groupid']

if rule in cluster_json:
    for key in cluster_json[rule]:
        params[key] = cluster_json[rule][key]
elif rule in custom_config_rules:
    # specifify custom configureations specific rules
    if rule == 'calculate_cluster':
        if job_properties['wildcards']['partition'] == 'onlyWELL':
            params['partition']='quick'
            params['time'] = '2:00:00'
            params['mem'] = '50G'
        elif job_properties['wildcards']['method'] == 'CCA':
            params['mem'] = '300G'
            params['partition']='norm'
            params['time'] = '12:00:00'
        else:
            params['partition']='norm'
            params['mem'] = '200G'
    if rule == 'calculate_umap':
        if job_properties['wildcards']['partition'] == 'onlyWELL':
            params['partition']='quick'
            params['time'] = '2:00:00'
            params['mem'] = '50G'
        else:
            params['time'] = '8:00:00'
            params['mem'] = '200G'
            params['partition']='norm'
            params['mem'] = '200G'
    if rule == 'perf_metrics':
        if job_properties['wildcards']['partition'] == 'onlyWELL':
            params['mem'] = '40G'
            params['partition']='quick'
            params['time'] = '1:00:00'
        elif job_properties['wildcards']['method'] == 'CCA':
            params['mem'] = '400G'
            params['partition']='largemem'
            params['time'] = '12:00:00'
        else:
            params['mem'] = '750G'
            params['partition']='largemem'
            params['time'] = '4:00:00'
    if rule == 'integrate_00':
        if re.search('scVI', job_properties['wildcards']['method']):
            params['partition']='gpu'
            params['extra'] = '--gres=gpu:1,lscratch:5'
            params['time'] = '8:00:00'
            params['mem'] = '100G'
        if job_properties['wildcards']['method'] == 'scVI':
            params['partition']='gpu'
            params['extra'] = '--gres=gpu:1,lscratch:5'
            params['time'] = '8:00:00'
            params['mem'] = '100G'
        elif job_properties['wildcards']['method'] == 'scVIprojection':
            params['partition']='gpu'
            params['extra'] = '--gres=gpu:k80:1,lscratch:20'
            params['time'] = '8:00:00'
            params['mem'] = '100G'
        elif job_properties['wildcards']['method'] == 'ldvae':
            params['partition']='gpu'
            params['extra'] = '--gres=gpu:1,lscratch:5'
            params['time'] = '48:00:00'
            params['mem'] = '200G'
        elif job_properties['wildcards']['method'] == 'CCA':
            params['partition']='largemem'
            params['time']='96:00:00'
            params['mem']='600G'
            params['time']='96:00:00'
            params['extra'] = '--gres=lscratch:5'
        elif job_properties['wildcards']['method'] == 'CCA' and job_properties['wildcards']['partition'] == 'onlyWELL':
            params['partition']='norm'
            params['time']='24:00:00'
            params['mem']= '50G'
            params['extra'] = '--gres=lscratch:1'
        elif job_properties['wildcards']['partition'] == 'onlyWELL':
            params['time']='2:00:00'
            params['mem']= '50G'
            params['extra'] = '--gres=lscratch:1'
        elif job_properties['wildcards']['transform'] == 'scran':
            params['mem']= '500G'
            params['partition']='largemem'
            params['time']='24:00:00'
            params['extra'] = '--gres=lscratch:5'
        else:
            print(job_properties['wildcards']['method'])
            params['partition']='gpu'
            params['extra'] = '--gres=gpu:1,lscratch:5'
            params['time'] = '8:00:00'
            params['mem'] = '200G'
    if rule == "make_seurat_objs" and job_properties['wildcards']['transform'] == 'scran':
        params['mem'] = '1000G'

else:# use default parameters
    params = cluster_json['__default__'] 

#%%
#allow for rules that have no wildcards 
if 'wildcards' in job_properties:
    ec_strings = [f"{key}={job_properties['wildcards'][key]}" for key in job_properties['wildcards'] ]
    ec_strings = '-'.join(ec_strings)
    output = f'{rule}.{ec_strings}'
else:
    #inlcude jobid for better error parsing
    jid=job_properties['jobid']
    output = f'{rule}-{jid}'


    
sbcmd=f'''sbatch --cpus-per-task={params['cpus-per-task']} \
    --mem={params['mem']} \
    --time={params['time']} \
    --job-name={rule} \
    --partition={params['partition']} \
    --output=00log/{output}.out \
    --error=00log/{output}.err \
    {params['extra']} \
    {jobscript}

'''
print(sbcmd)
os.system(sbcmd)

# %%
