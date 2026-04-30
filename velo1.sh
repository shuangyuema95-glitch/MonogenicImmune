#!/bin/bash
#SBATCH -J RNAvelo_case
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -n 4
#SBATCH -e error
#SBATCH -o out

RMSK=/public/home/yuxiaomingroup/msy/biosoft/annovar/humandb/hg38_rmsk.gtf
sample1=/public/home/yuxiaomingroup/msy/work/AID/juzhen/case
CR_gtf=/public/home/yuxiaomingroup/msy/biosoft/cellranger/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf
out=/public/home/yuxiaomingroup/msy/work/AID/juzhen/VeloOut
velocyto run10x -m $RMSK $sample1 $CR_gtf
