#!/bin/sh

#PBS -N per_batch
#PBS -k oe
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -q med-bio
i=${PBS_ARRAY_INDEX}

source ${PBS_O_WORKDIR}/aw_params.sh

module load plink R gcc

cd $TMPDIR

##############
### PRE-QC ###
##############

### preprocessing
# make BED
if [ "$ZCALLINPUT" = true ];then
  cp $(echo ${INFILES}|cut -d' ' -f${i}) n${i}.raw
  python ${PBS_O_WORKDIR}/convertReportToTPED.py -R n${i}.raw -O n${i}.tmp
  plink --tfile n${i}.tmp --make-bed --out n${i}.raw
else
  BASENAME=$(echo ${INFILES} | cut -d' ' -f${i} | sed 's/.bed$//g')
  cp ${BASENAME}.bed n${i}.raw.bed
  cp ${BASENAME}.bim n${i}.raw.bim
  cp ${BASENAME}.fam n${i}.raw.fam
fi

# update alleles and genome build
if [[ ! -z "$AFFYIDS" ]];then
  mv n${i}.raw.bed n${i}.bed
  mv n${i}.raw.bim n${i}.bim
  mv n${i}.raw.fam n${i}.fam
else
  plink --noweb --bfile n${i}.raw --update-alleles ${STRANDPATH}.update_alleles.txt --make-bed --out n${i}
  ${PBS_O_WORKDIR}/update_build.sh n${i} ${STRANDPATH}-${GENOMEBUILD}.strand n${i}
fi

### update sex with detected one
plink --bfile n${i} --check-sex --out n${i}
plink --bfile n${i} --update-sex n${i}.sexcheck 2 --make-bed --out n${i}

##############
### QC #######
##############

### individuals with excessive missing data and abnormal heterozygosity rate
plink --bfile n${i} --missing --out n${i}
plink --bfile n${i} --het --out n${i}
Rscript ${PBS_O_WORKDIR}/aw_imiss-vs-het.R n${i} ${HET} ${IMISS}
mv fail-imisshet-qc.txt n${i}.fail-imisshet-qc.txt
mv imiss-vs-het.png n${i}.imiss-vs-het.png
plink --bfile n${i} --remove n${i}.fail-imisshet-qc.txt --make-bed --out n${i}

cp n${i}.{bed,bim,fam} ${OUTPATH}/batch_data/
cp n${i}.{fail-imisshet-qc.txt,imiss-vs-het.png} ${OUTPATH}/report/data/
