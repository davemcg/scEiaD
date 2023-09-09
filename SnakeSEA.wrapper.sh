#!/bin/bash

# to run snakemake as batch job
# module load snakemake || exit 1

# activate conda
source /data/mcgaugheyd/conda/etc/profile.d/conda.sh

mkdir -p 00log 
sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"

# Run once to generate "well_batches.txt" and "droplet_counts.txt" for the SnakeSEA
# yes this is janky, but making it "clean" is not worth my time right now
#/data/mcgaugheyd/conda/envs/seacells/bin/python /home/mcgaugheyd/git/scEiaD/src/make_seacells.pull.py

snakemake -s /home/mcgaugheyd/git/scEiaD/SnakeSEA \
-pr --jobs 1999 \
--configfile $1 \
--use-conda \
--cluster-config /home/mcgaugheyd/git/scEiaD/cluster.json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 0 \
--resources parallel=4 


conda deactivate
