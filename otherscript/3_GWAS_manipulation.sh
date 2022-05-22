#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-330
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

# Fix allele oritation, remove error alleles for EA, and calculate meta GWAS
idx=$SLURM_ARRAY_TASK_ID

source /home1/zeyunlu/init.sh

echo Starting script 5 - fix allele orientation

params=`sed "${idx}q;d" /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/param.tsv`
echo "${idx} ${params}"
set -- junk $params
shift

CHR=$1
PHEN=$4

# make meta GWAS
Rscript /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/meta_gwas.R ${PHEN} ${CHR}

# make sig regions
Rscript /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/sig_regions.R ${PHEN} ${CHR}

echo DONE!
