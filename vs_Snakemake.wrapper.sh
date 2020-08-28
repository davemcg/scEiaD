#!/bin/bash

# to run snakemake as batch job
# module load snakemake || exit 1

# activate conda
mkdir -p 00log 
snakefile=$1
config=$2
cluster_json$3
cluster_config=$4
conda_sh=$5

source ${conda_sh}

snakemake -s $snakefile \
-pr --jobs 1999 \
--configfile $config \
--use-conda \
--cluster "./${cluster_config} ${cluster_json}"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 0 \
--resources parallel=4

