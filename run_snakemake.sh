#!/bin/bash
PROJECT="$1" snakemake --cluster-config clusterconfig.yaml \
--cluster "qsub -S /bin/bash -R y -j y -cwd -V -m a -M zzj.zju@gmail.com -l h_data={cluster.h_data},h_rt={cluster.h_rt} -pe shared {cluster.threads}" \
--latency-wait 60 --jobs 100 \
--rerun-incomplete \
$2
