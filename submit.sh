#!/bin/bash
#$ -S /bin/bash
#$ -l h_data=2G,h_rt=24:00:0
#$ -R y
#$ -V
#$ -cwd
#$ -j y
#$ -m be
#$ -M zzj.zju@gmail.com
#


source activate py36
./run_snakemake.sh $1
