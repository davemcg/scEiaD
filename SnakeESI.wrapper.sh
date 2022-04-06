#!/bin/bash

# to run snakemake as batch job
# module load snakemake || exit 1

# activate conda
#rm -rf 00log/
mkdir -p 00log/ 
snakefile=$1
config=$2
cluster_json=$3
cluster_config=$4
conda_sh=$5
# activate conda
source ${conda_sh}
module load python
snakemake -s $snakefile \
-pr --jobs 1999 \
--configfile $config \
--use-conda \
--cluster "python3 ${cluster_config} ${cluster_json}"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 0 \
--resources parallel=4 \
--resources integration_gpu=8

