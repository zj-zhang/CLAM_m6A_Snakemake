#!/bin/bash
PROJECT_NAME="$1" snakemake --cluster-config config.yaml \
--cluster "qsub -R y -j y -cwd -V -m a -M zzj.zju@gmail.com -l h_data={cluster.h_data},h_rt={cluster.h_rt}" \
--latency-wait 60 --jobs 4 \
--rerun-incomplete \
$2
