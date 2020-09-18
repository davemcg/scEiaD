#!/bin/bash

# to run snakemake as batch job
# module load snakemake || exit 1

# activate conda
source /data/mcgaugheyd/conda/etc/profile.d/conda.sh
conda activate kallisto_046

mkdir -p 00log 
sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s /home/mcgaugheyd/git/massive_integrated_eye_scRNA/SnakePOP \
-pr --jobs 1999 \
--configfile /home/mcgaugheyd/git/massive_integrated_eye_scRNA/config.yaml \
--use-conda \
--cluster-config /home/mcgaugheyd/git/massive_integrated_eye_scRNA/cluster.json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 0 \
--resources parallel=4


conda deactivate
