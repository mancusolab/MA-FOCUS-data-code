#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-3400

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

DATAF=/project/nmancuso_8/zeyunlu/data/mefocus/sim_genos/eur_afr_eas
SCRIPTF=/project/nmancuso_8/zeyunlu/projects/mefocus/scripts

start_time=$SECONDS
params=`sed "${NR}q;d" ${SCRIPTF}/param_eur_afr_eas.tsv`
echo "${NR} ${params}"
set -- junk $params
shift

# PARAMETERS
ROW=$1
LOCUS=$2
POP=$3
MODEL=$4
eQTLMODEL=$5
H2G=$6
H2GE=$7
N1=$8
NGE1=$9
N2=${10}
NGE2=${11}
N3=${12}
NGE3=${13}

TMP=/scratch/zeyunlu/tmp/mefocus/eur_afr_eas/

echo "Running on Locus ${LOCUS} on Row ${ROW}"

OUT=$TMP/row${ROW}.locus${LOCUS}

python ${SCRIPTF}/mefocus_run_eur_afr_eas.py \
$DATAF/locus${LOCUS}/${LOCUS}.EUR \
$DATAF/locus${LOCUS}/${LOCUS}.AFR \
$DATAF/locus${LOCUS}/${LOCUS}.EAS \
--ngwas_pop1 $N1 \
--nqtl_pop1 $NGE1 \
--ngwas_pop2 $N2 \
--nqtl_pop2 $NGE2 \
--ngwas_pop3 $N3 \
--nqtl_pop3 $NGE3 \
--model $MODEL \
--share-model $eQTLMODEL \
--eqtl-h2 $H2G \
--var-explained $H2GE \
--output $OUT \
--sim $ROW \
--locus $LOCUS \
--seed `python -c "print( int(int($LOCUS) + 1124))"` \
--pop $POP

if [[ $LOCUS != 1 ]]
then
  tail -n+2 ${OUT}.focus.tsv > ${OUT}.focus.nocol.tsv
  rm -rf ${OUT}.focus.tsv
fi
