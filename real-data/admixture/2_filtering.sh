#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=8Gb

# Filter out indels, missing data, rare variants
dir="/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/"
in_data="1000G_GENOA_snpsOnly_rawMerge"
module load gcc/8.3.0 
module load plink/1.07

# Filter for missingness and minor allele frequency
cd $dir
#plink2 --bfile $in_data --geno 0.01 --maf 0.01 --max-maf 0.99 --make-bed --out 1000G_GENOA_snpsOnly_geno0.01_maf0.01

# Update population labels
#cat /project/nmancuso_8/data/GENOA/raw/genotype/aa/rename_AA.list > update_pop_IDs.list
#cat /project/nmancuso_8/data/GENOA/raw/genotype/ea/rename_EA.list >> update_pop_IDs.list

#fam="1000G_GENOA_snpsOnly_geno0.01_maf0.01.fam"
#key="/project/nmancuso_8/sgopalan/scripts/GENOA_ADMIXTURE/misc_files/igsr_samples.tsv"
#N=`wc -l $fam | cut -f1 -d' '`
#for i in `seq 1 $N`;do
#  ID=`sed -n ${i}p $fam | cut -f2`
#  POP=`grep -w $ID $key | cut -f4`
#  if [ "$POP" != "" ];then
#    echo 0 $ID $POP $ID >> update_pop_IDs.list
#  fi
#done
#plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01 --update-ids update_pop_IDs.list --make-bed --out 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated

# Find SNPs out of HWE and make list of unrelateds
POPS=`awk '{print $1}' 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated.fam | sort | uniq`
#for POP in ${POPS[@]};do
#  grep $POP 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated.fam > ${POP}.list
#  plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated --keep ${POP}.list --hardy midp --out $POP
#  awk '{if ($10<0.00001) {print $2}}' ${POP}.hardy >> HWE0.00001_varsToRemove.list
#  rm -f ${POP}.*
#done
#plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated --exclude HWE0.00001_varsToRemove.list --make-bed --out 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001

# Remove all AT/CG SNPs
#awk '{if ($5=="A" && $6=="T" || $5=="T" && $6=="A" || $5=="C" && $6=="G" || $5=="G" && $6=="C"){print $2}}' 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001.bim > AT-CG_SNPs.list
#plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001 --exclude AT-CG_SNPs.list --make-bed --out 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG

# LD filtering
#plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG --indep-pairwise 200kb 0.3 --out LD_200kb-1kbwin-0.3
#plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG --exclude LD_200kb-1kbwin-0.3.prune.out --make-bed --out 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3

# Clean up
if [ -f 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3.bed ];then
  rm -f 1000G_GENOA_snpsOnly_geno0.01_maf0.01.*
  rm -f 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated.*
  rm -f 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001.*
  rm -f 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG.*
fi
