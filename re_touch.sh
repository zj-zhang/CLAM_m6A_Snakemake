#!/bin/bash
# re-touch files by generating order
# so to avoid snakemake re-run from start
# Zijun Zhang
# 1.3.2018
project=$1
touch $project/config/*
find $project/star/ -exec touch {} \;
find $project/macs2/ -exec touch {} \;
touch $project/clam/peaks*/*unique*
touch $project/clam/peaks*/*combine*
touch $project/clam/peaks*/*rescue*
touch $project/clam/peaks*/*macs2*
touch $project/clam/peaks*/*.png
find $project/homer/ -exec touch {} \;
find $project/repeats/ -exec touch {} \;
find $project/topology/ -exec touch {} \;
find $project/reports/ -exec touch {} \;
find $project/bigwig/ -exec touch {} \;
find $project/archive/ -exec touch {} \;
