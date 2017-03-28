module load plink

PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/genotype_QC_lite/aw_params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/

cd $CMDPATH
qsub ${CMDPATH}/aw_step2.sh


