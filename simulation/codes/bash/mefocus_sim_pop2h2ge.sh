#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-1400

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

DATAF=/project/nmancuso_8/zeyunlu/data/mefocus/sim_genos/eur_afr
SCRIPTF=/project/nmancuso_8/zeyunlu/projects/mefocus/scripts

start_time=$SECONDS
params=`sed "${NR}q;d" ${SCRIPTF}/param_pop2h2ge.tsv`
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
N1=$7
NGE1=$8
N2=$9
NGE2=${10}
H2GE1=${11}
H2GE2=${12}

TMP=/scratch/zeyunlu/tmp/mefocus/pop2h2ge/

echo "Running on Locus ${LOCUS} on Row ${ROW}"

OUT=$TMP/row${ROW}.locus${LOCUS}

python ${SCRIPTF}/mefocus_run_pop2h2ge.py \
$DATAF/locus${LOCUS}/${LOCUS}.EUR \
$DATAF/locus${LOCUS}/${LOCUS}.AFR \
--ngwas_pop1 $N1 \
--nqtl_pop1 $NGE1 \
--ngwas_pop2 $N2 \
--nqtl_pop2 $NGE2 \
--model $MODEL \
--share-model $eQTLMODEL \
--eqtl-h2 $H2G \
--var-explained_p1 $H2GE1 \
--var-explained_p2 $H2GE2 \
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
