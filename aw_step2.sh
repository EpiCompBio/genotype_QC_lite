#!/bin/sh

#PBS -N zcall
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

cd $TMPDIR
cp ${ZCALLPATH}/n${i}.zcall .

### preprocessing
python ${PBS_O_WORKDIR}/convertReportToTPED.py -R n${i}.zcall -O n${i}.tmp
plink --tfile n${i}.tmp --make-bed --out n${i}.raw
plink --noweb --bfile n${i}.raw --update-alleles ${STRANDPATH}.update_alleles.txt --make-bed --out n${i}
${PBS_O_WORKDIR}/update_build.sh n${i} ${STRANDPATH}-${GENOMEBUILD}.strand n${i}

### sex problems
plink --bfile n${i} --check-sex --out n${i}
grep PROBLEM n${i}.sexcheck > n${i}.sexprobs
# TODO write to fail-sexcheck-qc.txt

### missing data and abnormal heterozygosity rate
plink --bfile n${i} --missing --out n${i}
plink --bfile n${i} --het --out n${i}
# TODO check that all required R packages are installed, including geneplotter
Rscript ${PBS_O_WORKDIR}/imiss-vs-het_modified.R n${i}

### relatedness
# TODO do this once batches are merged?
plink --bfile n${i} --indep-pairwise 50 5 0.2 --out n${i}
plink --bfile n${i} --extract n${i}.prune.in --genome --out n${i}
Rscript ${PBS_O_WORKDIR}/plot-IBD_modified.R n${i}
# TODO write to fail-IBD-qc.txt

find . -type f |grep -vE '(.zcall$|.tmp.|.py$)' |xargs -ILIST cp LIST ${PLINKPATH}
