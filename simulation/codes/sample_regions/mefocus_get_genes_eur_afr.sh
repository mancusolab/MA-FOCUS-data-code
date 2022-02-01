#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-300

# tdir is where the intermediate plink files are put
if [ ! $SLURM_ARRAY_TASK_ID ]; then
  idx=$1
  tdir=./
else
  idx=$SLURM_ARRAY_TASK_ID
  tdir=$SCRATCHDIR
fi

plink=/project/nmancuso_8/zeyunlu/tools/plink2

# gene complexity of sampled region
lower=5
upper=8

DATAF=/project/nmancuso_8/zeyunlu/data/mefocus/
SCRIPTF=/project/nmancuso_8/zeyunlu/projects/mefocus/scripts/

# set up variables for reference data needed
hapmap=${DATAF}HAPMAP_SNPS/
ind_loci=${DATAF}EUR_AFR_LDETECT_indLoci.bed
genes=${DATAF}glist-hg19.nodupe.autosome

RESF=${DATAF}/sim_genos/eur_afr/
if [ ! -d ${RESF} ]
then
  mkdir -p ${RESF}
fi

odir=${RESF}/locus$idx
if [ -d $odir ]; then
  rm -rf $odir
fi
mkdir -p $odir

python ${SCRIPTF}sample_genes.py $ind_loci $genes \
-l $lower -u $upper -o ${odir}/gene.${idx}.list \
--loc_output ${odir}/locus.${idx}.list \
--seed `python -c "print( int(int($idx) + 1124))"`

params=`sed "1q;d" ${odir}/locus.${idx}.list`
set -- junk $params
shift

# get the region bounds
chr=$1
numchr=`echo $chr | sed 's/chr//'`
locus_start=$2
locus_stop=$3

for megroup in EUR AFR
do
  # replace with your path to PLINK reference data
  $plink --bfile ${DATAF}/1000GP3/${megroup}/1000G.${megroup}.QC.${numchr} --chr $numchr --from-bp $locus_start --to-bp $locus_stop --make-bed \
  --out $odir/${idx}.${megroup} --snps-only --hwe midp 1e-5 --geno 0.01 --maf 0.01 \
  --memory 2048 --extract $hapmap/hm.$numchr.snp --silent --force-intersect
done
