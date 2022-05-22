#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-15
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

# filter out large chisq value and split by chr

idx=$SLURM_ARRAY_TASK_ID

source /home1/zeyunlu/init.sh

params=`sed "${idx}q;d" /home1/zeyunlu/research/sub_mefocus/scripts/param3.tsv`
  set -- junk $params
  shift
PHEN=$1

Rscript /home1/zeyunlu/research/sub_mefocus/scripts/filter_chisq.R ${PHEN}
