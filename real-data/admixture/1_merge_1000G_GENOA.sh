#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH --mem=128Gb
#SBATCH --partition=oneweek

source ~/.bashrc

# Merge 1000G autosomal data
outdir="/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/"
scripts_dir="/project/nmancuso_8/sgopalan/scripts/GENOA_ADMIXTURE/"
#if [ ! -f ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz ];then
#  for chr in `seq 1 22`;do
#    file=/project/nmancuso_8/data/1000GP3_RESEQ/vcf/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
#    if [ $chr -eq 1 ];then
#      echo $file > ${scripts_dir}misc_files/merge_chr.list
#    else
#      echo $file >> ${scripts_dir}misc_files/merge_chr.list
#    fi
#  done
#  module load gcc/8.3.0
#  module load bcftools/1.9
#  bcftools concat -f ${scripts_dir}misc_files/merge_chr.list > ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf
#  bgzip -f ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf
#fi

# Recode vcf version
#zcat ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | head -n1 > ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf
#sed -i "0,/4.3/{s/4.3/4.2/;}" ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf
#zcat ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | tail -n +2 >> ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf
#bgzip ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf
#if [ -f ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf.gz ];then
#  rm -f ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
#fi
#tabix ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf.gz

# Extract on certain individuals from 1000G (African, European, and American populations)
IDs_file="${scripts_dir}misc_files/igsr_samples.tsv"
data=${outdir}1000G_AFR-EUR-AMR-pops_allvariants.vcf.gz
#if [ ! -f ${data}.tbi ];then
#  module load gcc/8.3.0
#  module load vcftools
#  vcftools --gzvcf ${outdir}ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_recoded.vcf.gz --keep $IDs_file --recode --stdout | bgzip -c > ${data}
#  tabix $data
#fi

# Convert to plink binary
GENOA_data="/project/nmancuso_8/data/GENOA/processed/genotype/gwas_datasets/ea-aa_intersectMerge_R2-0.6_MAF0.01_HWE0.000001_filtered_no-wonky-vars.vcf"
module load gcc/8.3.0
module load plink/1.07
#plink2 --vcf $data --make-bed --out ${outdir}1000G_AFR-EUR-AMR-pops_allvariants
#plink2 --vcf $GENOA_data --make-bed --out ${outdir}GENOA_data

# Rename variants
#for file in 1000G_AFR-EUR-AMR-pops_allvariants GENOA_data;do
#  plink1.9 --bfile ${outdir}${file} --set-missing-var-ids  @:'#'\$1,\$2 --snps-only --make-bed --out ${outdir}${file}_snps-renamed
#done

# Merge 1000G and GENOA data
#plink1.9 --bfile ${outdir}GENOA_data_snps-renamed --bmerge ${outdir}1000G_AFR-EUR-AMR-pops_allvariants_snps-renamed --make-bed --out ${outdir}1000G_GENOA_snpsOnly_rawMerge
#plink1.9 --bfile ${outdir}GENOA_data_snps-renamed --exclude ${outdir}1000G_GENOA_snpsOnly_rawMerge-merge.missnp --make-bed --out ${outdir}GENOA_data_snps-renamed_clean
#plink1.9 --bfile ${outdir}1000G_AFR-EUR-AMR-pops_allvariants_snps-renamed --exclude ${outdir}1000G_GENOA_snpsOnly_rawMerge-merge.missnp --make-bed --out ${outdir}1000G_AFR-EUR-AMR-pops_allvariants_snps-renamed_clean
plink1.9 --bfile ${outdir}GENOA_data_snps-renamed_clean --bmerge ${outdir}1000G_AFR-EUR-AMR-pops_allvariants_snps-renamed_clean --make-bed --out ${outdir}1000G_GENOA_snpsOnly_rawMerge
