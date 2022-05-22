#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1-30
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

idx=$SLURM_ARRAY_TASK_ID

source /home1/zeyunlu/init.sh

params=`sed "${idx}q;d" /home1/zeyunlu/research/sub_mefocus/scripts/param2.tsv`
  set -- junk $params
  shift
PHEN=$1
POP=$2

cd /project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/${POP}_${PHEN}_GWAMA
# rm -rf ${POP}_${PHEN}_sumstats_new.*
#
# for CHR in `seq 1 22`
# do
#     if [ $CHR -eq 1 ]
#     then
#       zcat ${POP}_${PHEN}_chr${CHR}_GWAMA_GRCh38_dbSNP-updated_wHeader.tsv.gz | cut -f1-12,14- > ${POP}_${PHEN}_new.tsv
#     else
#       zcat ${POP}_${PHEN}_chr${CHR}_GWAMA_GRCh38_dbSNP-updated_wHeader.tsv.gz | tail -n+2 | cut -f1-12,14- >> ${POP}_${PHEN}_new.tsv
#     fi
#   done
#
# gzip -f ${POP}_${PHEN}_new.tsv

source activate ldsc

python2 /project/nmancuso_8/zeyunlu/tools/ldsc/munge_sumstats.py \
  --sumstats ${POP}_${PHEN}_new.tsv.gz \
  --N-col n_samples \
  --snp rs_number \
  --a1 reference_allele \
  --a2 other_allele \
  --p p-value \
  --frq eaf \
  --signed-sumstats z,0 \
  --out ${POP}_${PHEN}_new_munged
