#!/bin/sh

module load plink
module load R
# TODO check it's installed: devtools::install_github("gabraham/flashpca/flashpcaR")
# TODO check that all required R packages are installed, including geneplotter

PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/genotype_QC_lite/aw_params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/
mkdir -p ${PLINKPATH}

cd ${CMDPATH}
qsub aw_step2.sh
#qsub aw_step3.sh

