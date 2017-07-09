#!/bin/sh

### load parameters
PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/genotype_QC_lite/aw_params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/
mkdir -p ${OUTPATH}/report/data
mkdir -p ${OUTPATH}/batch_data

### load software and packages
module load plink R gcc
Rscript ${CMDPATH}/aw_check_packages.r

### run analysis
cd ${CMDPATH}
NSAMPLES=$(echo $INFILES | wc -w)
if [ $NSAMPLES == "1" ]; then NSAMPLES=2:2; fi
per_batch=$(qsub -J 1-${NSAMPLES} aw_per_batch.sh)
all_batches=$(qsub -W depend=afterok:$per_batch aw_all_batches.sh)

