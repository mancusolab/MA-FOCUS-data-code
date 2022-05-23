#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1-330
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

idx=$SLURM_ARRAY_TASK_ID

source /home1/zeyunlu/init.sh

params=`sed "${idx}q;d" /project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/param.tsv`
echo "${idx} ${params}"
set -- junk $params
shift

CHR=$1
PHEN=$4

fdir=/project/nmancuso_8/zeyunlu/projects/sub_mefocus/scripts/
wdir=/project/nmancuso_8/sgopalan/results/FUSION_GENOA/
lddir=/project/nmancuso_8/data/1000GP3_multiPop_allelesAligned/

# EA
ea_ss=/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/EA/EA_${PHEN}_new_munged_chr${CHR}.tsv.gz
ea_out=/project/nmancuso_8/data/MEFOCUS/twas_results/EA/${PHEN}_twas_ea_chr${CHR}.tsv

Rscript $fdir/twas_test.R \
--sumstats ${ea_ss} \
--out ${ea_out} \
--weights ${wdir}ea_wgt_pos_file.tsv \
--weights_dir $wdir/ea/ \
--ref_ld_chr $lddir/EUR/1000G.EUR.QC.allelesAligned. \
--chr $CHR \
--trait $PHEN

# AA
aa_ss=/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/AA/AA_${PHEN}_new_munged_chr${CHR}.tsv.gz
aa_out=/project/nmancuso_8/data/MEFOCUS/twas_results/AA/${PHEN}_twas_aa_chr${CHR}.tsv

Rscript $fdir/twas_test.R \
--sumstats ${aa_ss} \
--out ${aa_out} \
--weights ${wdir}aa_wgt_pos_file.tsv \
--weights_dir $wdir/aa/ \
--ref_ld_chr $lddir/AFR/1000G.AFR.QC.allelesAligned. \
--chr $CHR \
--trait $PHEN

# Meta
meta_ss=/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/meta/meta_${PHEN}_new_munged_chr${CHR}.tsv.gz
meta_out=/project/nmancuso_8/data/MEFOCUS/twas_results/meta/${PHEN}_twas_meta_chr${CHR}.tsv

Rscript $fdir/twas_test.R \
--sumstats ${meta_ss} \
--out ${meta_out} \
--weights ${wdir}ea_wgt_pos_file.tsv \
--weights_dir $wdir/ea/ \
--ref_ld_chr $lddir/EUR/1000G.EUR.QC.allelesAligned. \
--chr $CHR \
--trait $PHEN
