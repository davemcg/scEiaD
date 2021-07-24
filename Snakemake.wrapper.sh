#!/bin/bash

# run in the data folder for this project
# on biowulf2:
# /data/mcgaugheyd/projects/nei/brooks/oca_rna-seq

mkdir -p 00log

module load snakemake || exit 1

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s $1 \
  -pr --local-cores 2 --jobs 500 \
  --cluster-config $3 \
  --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
  --configfile $2 --use-conda \
  -k --restart-times 0 
