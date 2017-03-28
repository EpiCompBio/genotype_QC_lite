#!/bin/sh

#PBS -N per_batch
#PBS -k oe
#PBS -e ./logs/
#PBS -o ./logs/
#PBS -l mem=2gb
#PBS -l walltime=2:00:00
#PBS -q med-bio
#PBS -J 1-2
i=${PBS_ARRAY_INDEX}

source ${PBS_O_WORKDIR}/aw_params.sh

module load plink
module load R
# TODO check that all required R packages are installed, including geneplotter

cd $TMPDIR
cp ${ZCALLPATH}/n${i}.zcall .

##############
### PRE-QC ###
##############

### preprocessing
python ${PBS_O_WORKDIR}/convertReportToTPED.py -R n${i}.zcall -O n${i}.tmp
plink --tfile n${i}.tmp --make-bed --out n${i}.raw
plink --noweb --bfile n${i}.raw --update-alleles ${STRANDPATH}.update_alleles.txt --make-bed --out n${i}
${PBS_O_WORKDIR}/update_build.sh n${i} ${STRANDPATH}-${GENOMEBUILD}.strand n${i}

### update sex with detected one
plink --bfile n${i} --check-sex --out n${i}
plink --bfile n${i} --update-sex n${i}.sexcheck 2 --make-bed --out n${i}

##############
### QC #######
##############

### missing data and abnormal heterozygosity rate
plink --bfile n${i} --missing --out n${i}
plink --bfile n${i} --het --out n${i}
Rscript ${PBS_O_WORKDIR}/imiss-vs-het_modified.R n${i}
plink --bfile n${i} --remove n${i}.fail-imisshet-qc.txt --make-bed --out n${i}

### markers with excessive missing data
#plink --bfile n${i}.clean-inds --missing --out n${i}.clean-inds
#Rscript ../lmiss-hist_modified.R n1
plink --bfile n${i} --geno 0.05 --write-snplist --out n${i}

cp n${i}.{bed,bim,fam,snplist} ${PLINKPATH}
