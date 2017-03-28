module load plink

PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/

cd CMDPATH
qsub ${CMDPATH}/aw_step2.sh


