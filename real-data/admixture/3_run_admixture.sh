#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --array=1-120

source ~/.bashrc
i=$SLURM_ARRAY_TASK_ID

d="/project/nmancuso_8/sgopalan/results/GENOA_ADMIXTURE/"
r="admixture_results/"

# Set parameters
FILE="1000G_GENOA_snpsOnly_geno0.01_maf0.01_popUpdated_HWE0.00001_noAT-CG_LD-LD_200kb-1kbwin-0.3"
K=`sed -n ${i}p /project/nmancuso_8/sgopalan/scripts/GENOA_ADMIXTURE/misc_files/admixture.params | awk '{print $1}'`
rep=`sed -n ${i}p /project/nmancuso_8/sgopalan/scripts/GENOA_ADMIXTURE/misc_files/admixture.params | awk '{print $2}'`

# Create directory structure
cd $d
mkdir -p $r/admixture_rep${rep}

# Run admixture
cd $r/admixture_rep${rep}
echo rep K seed log-likelihood cross-val_error > run_stats.txt
admixture -s time --cv ${d}${FILE}.bed $K | tee log_rep${rep}_K${K}.out

# Rename output files
for ext in Q P;do
  mv ${FILE}.${K}.${ext} ${d}${r}K${K}_rep${rep}.${ext}
done

# Extract info from the run
SEED=`grep seed log_rep${rep}_K${K}.out|cut -f 3 -d ' '`
LLIKE=`grep ^Log log_rep${rep}_K${K}.out|cut -f 2 -d ' '`
CVE=`grep CV log_rep${rep}_K${K}.out|cut -f 4 -d ' '`
echo $rep $K $SEED $LLIKE $CVE >> run_stats.txt

echo DONE!
