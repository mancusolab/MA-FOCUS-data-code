#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=40Gb
#SBATCH --array=1-330
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

source /home1/zeyunlu/init.sh

idx=$SLURM_ARRAY_TASK_ID

params=`sed "${idx}q;d" /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/param.tsv`
echo "${idx} ${params}"
set -- junk $params
shift

CHR=$1
PHEN=$4
MINP=1

R=/project/nmancuso_8/zeyunlu/tools/miniconda3/bin/Rscript

LDETECT_REGS=/project/nmancuso_8/data/LDETECT/EUR_AFR_LDETECT_indLoci.GRCh38.split.bed
DATA_DIR=/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/
LD_DIR=/project/nmancuso_8/data/1000GP3_multiPop_allelesAligned/

EA_TWAS=/project/nmancuso_8/data/MEFOCUS/twas_results/modified/EA_${PHEN}_twas_modified_chr${CHR}.tsv
AA_TWAS=/project/nmancuso_8/data/MEFOCUS/twas_results/modified/AA_${PHEN}_twas_modified_chr${CHR}.tsv
META_TWAS=/project/nmancuso_8/data/MEFOCUS/twas_results/modified/meta_${PHEN}_twas_modified_chr${CHR}.tsv

SIG_DIR=/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_sig_region/chr/
SIG_REGION=${SIG_DIR}EA_${PHEN}_sig_region_chr${CHR}.bed:${SIG_DIR}AA_${PHEN}_sig_region_chr${CHR}.bed:${SIG_DIR}meta_${PHEN}_sig_region_chr${CHR}.bed

# general mefocus
echo "Running ME-FOCUS"
$R /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/ma_focus.R \
  --TWAS ${EA_TWAS}:${AA_TWAS} \
  --chr $CHR \
  --regions $LDETECT_REGS \
  --sig_regions $SIG_REGION \
  --ref_ld_chr ${LD_DIR}EUR/1000G.EUR.QC.allelesAligned.${CHR}:${LD_DIR}AFR/1000G.AFR.QC.allelesAligned.${CHR} \
  --minp.input ${MINP}:${MINP} \
  --perform_single_pop_FOCUS T \
  --phen ${PHEN} \
  --prior_prob BLOCKGENE \
  --gencode /project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/gencode.v26.GRCh38.genes.only.fp.tsv \
  --out /project/nmancuso_8/data/MEFOCUS/focus_results/mefocus38_2/${PHEN}_chr${CHR}_mefocus.tsv

# general meta mefocus
echo "Running meta FOCUS"
$R /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/ma_focus.R \
  --TWAS ${META_TWAS} \
  --chr $CHR \
  --regions $LDETECT_REGS \
  --sig_regions $SIG_REGION \
  --ref_ld_chr ${LD_DIR}EUR/1000G.EUR.QC.allelesAligned.${CHR} \
  --minp.input ${MINP} \
  --perform_single_pop_FOCUS T \
  --phen ${PHEN} \
  --prior_prob BLOCKGENE \
  --gencode /project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/gencode.v26.GRCh38.genes.only.fp.tsv \
  --out /project/nmancuso_8/data/MEFOCUS/focus_results/mefocus38_2/${PHEN}_chr${CHR}_meta.tsv

echo "Finish Running"
