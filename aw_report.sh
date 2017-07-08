#!/bin/sh

#PBS -N report
#PBS -k oe
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=1:00:00
#PBS -q med-bio

PARAMS=${PBS_O_WORKDIR}/aw_params.sh
source $PARAMS

module load R gcc

## generate report
cd ${PBS_O_WORKDIR}
cp -r assets ${OUTPATH}/report/
Rscript -e "library(knitr); knit('report.Rhtml', output='${OUTPATH}/report/report.html')" $PARAMS

### housekeeping
chmod -R 750 $OUTPATH

