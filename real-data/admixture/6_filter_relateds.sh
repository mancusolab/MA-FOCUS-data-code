#!/bin/bash
cd /project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE

cut -f1,2 -d, AA_GENOA_0.75-AFR-anc.csv | tail -n+2 | sed "s/,/\t/g" > high-AFR-anc_AA_GENOA.list

# Kinship analyses
module load gcta/1.93
mkdir -p kinship_results
cd kinship_results
plink2 --bfile 1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3 --keep high-AFR-anc_AA_GENOA.list --make-bed --out high-AFR-anc_AA_GENOA
gcta64 --bfile high-AFR-anc_AA_GENOA --thread-num 48 --make-grm --out high-AFR-anc_AA_GENOA
gcta64 --grm high-AFR-anc_AA_GENOA --grm-cutoff 0.05 --thread-num 48 --make-grm --out high-AFR-anc_AA_GENOA_grm0.05_unrelateds
rm -f high-AFR-anc_AA_GENOA.bed high-AFR-anc_AA_GENOA.bim high-AFR-anc_AA_GENOA.fam
