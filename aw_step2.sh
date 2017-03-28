#!/bin/sh

#PBS -N zcall
#PBS -k oe
#PBS -e ./logs/
#PBS -o ./logs/
#PBS -l mem=2gb
#PBS -l walltime=2:00:00
#PBS -q med-bio
#PBS -J 1-52
i=${PBS_ARRAY_INDEX}

source {PBS_O_WORKDIR}/aw_params.sh

cd $TMPDIR
cp {PBS_O_WORKDIR}/convertReportToTPED.py ${ZCALLPATH}/n${i}.zcall .

python convertReportToTPED.py -R n${i}.zcall -O n${i}.tmp

module load plink
plink --tfile n${i}.tmp --make-bed --out n${i}.raw
plink --noweb --bfile n${i}.raw --update-alleles ${STRANDPATH}.update_alleles.txt --make-bed --out n${i}
${PBS_O_WORKDIR}/update_build.sh n${i} ${STRANDPATH}-${GENOMEBUILD}.strand n${i}

cp n${i}.* ${PLINKPATH}
