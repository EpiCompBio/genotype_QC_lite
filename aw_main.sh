#!/bin/sh

### load parameters
PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/genotype_QC_lite/aw_params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/
mkdir -p ${OUTPATH}

### load software and packages
module load plink R gcc
Rscript ${CMDPATH}/aw_check_packages.r

### run analysis
cd ${CMDPATH}
qsub aw_per_batch.sh
#qsub aw_all_batches.sh

