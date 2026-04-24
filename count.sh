#!/bin/bash
#SBATCH -J case_cellranger
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -n 4
#SBATCH -e error
#SBATCH -o out

ref=/public/home/yuxiaomingroup/msy/biosoft/cellranger/ref/refdata-gex-GRCh38-2020-A 
cr=/public/home/yuxiaomingroup/msy/biosoft/cellranger8/cellranger-8.0.1/cellranger # convenient for RNA velocity analysis
fo1=$PWD

samples=$(ls *.gz | cut -d"_" -f1 | sort -u | awk '{printf "--sample=%s ", $1}')

$cr count --id=case \
   --fastqs=$fo1 \
   $samples \
   --transcriptome=$ref \
   --create-bam true \
   --nosecondary

