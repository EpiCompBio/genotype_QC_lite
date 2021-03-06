#!/bin/sh

#PBS -N all_batches
#PBS -k oe
#PBS -l select=1:ncpus=10:mem=31gb
#PBS -l walltime=8:00:00
#PBS -q med-bio

PARAMS=${PBS_O_WORKDIR}/aw_params.sh
source $PARAMS

module load plink R gcc
PLINK="plink --threads 10"

######################
### MERGE ############
######################

### create single dataset
# only samples that pass missingness and heterogeneity filters
cd ${OUTPATH}/batch_data
ls | egrep '^n[0-9]+.bed$' | sed 's/.bed//g' | sort -n -k1.2 > ${TMPDIR}/batches.list
$PLINK --merge-list ${TMPDIR}/batches.list --make-bed --out ${TMPDIR}/all --merge-mode 2
cd ${TMPDIR}

### if needed, map Affy IDs to RS IDs
if [[ ! -z "$AFFYIDS" ]];then
  sed -e 's/\"//g' ${AFFYIDS} > mapping.txt
  $PLINK --bfile all --update-name mapping.txt 3 1 20 --make-bed --out all
fi

######################
### QC INDIVIDUALS ###
######################

### create smaller, less redundant dataset
$PLINK --bfile all --exclude ${HIGHLDFILE} --range --indep-pairwise 50 5 0.2 --out all

### IBD individuals
$PLINK --bfile all --extract all.prune.in --geno ${GENO} --hwe ${HWE} --out all --rel-cutoff ${REL}
cut all.fam -f1-2 -d' '|sort|uniq|sed -e 's/ /\t/g' > temp1
sort all.rel.id|uniq > temp2
comm -23 temp1 temp2 > fail-ibd-qc.txt

### update HapMap data to the same genome build
mkdir -p hapmap
HAPMAPBFILENEW=hapmap/$(basename $HAPMAPBFILE).${GENOMEBUILD}
$PLINK --bfile ${HAPMAPBFILE} --update-alleles ${STRANDPATH}.update_alleles.txt \
      --make-bed --out ${HAPMAPBFILENEW}
${PBS_O_WORKDIR}/update_build.sh ${HAPMAPBFILE} ${STRANDPATH}-${GENOMEBUILD}.strand ${HAPMAPBFILENEW}

### merge with HapMap data
cut -f2 all.bim | sort > all.sort
cut -f2 ${HAPMAPBFILENEW}.bim | sort > ${HAPMAPBFILENEW}.sort
comm -12 all.sort ${HAPMAPBFILENEW}.sort > hapmap_common.snps
$PLINK --bfile all --extract hapmap_common.snps --make-bed --out all.hapmap_snps
$PLINK --bfile all.hapmap_snps --bmerge ${HAPMAPBFILENEW} --out all.missnp
$PLINK --bfile ${HAPMAPBFILENEW} --flip all.missnp.missnp --make-bed --out ${HAPMAPBFILENEW}.flipped
#cat fail-ibd-qc.txt ${OUTPATH}/n*.fail-* | sort -k1 | uniq | cut -f1,2 > fail-qc-inds_preancestry.txt
$PLINK --bfile all.hapmap_snps --bmerge ${HAPMAPBFILENEW}.flipped --extract all.prune.in \
      --make-bed --out all.shared_hapmap_pruned

### individuals with divergent ancestry
Rscript ${PBS_O_WORKDIR}/aw_ancestry.r ${ANC} ${ETHNICFILE}

### remove individuals not passing dataset-wide QC
cat fail-* ${OUTPATH}/report/data/n*.fail-* | sort -k1 | uniq | cut -f1,2 > fail-qc-inds.txt
$PLINK --bfile all --remove fail-qc-inds.txt --make-bed --out all.clean-inds

######################
### QC MARKERS #######
######################

### check missing rate
$PLINK --bfile all.clean-inds -missing --out all.clean-inds
Rscript ${PBS_O_WORKDIR}/aw_lmiss-hist.R all.clean-inds

### remove markers not passing dataset-wide QC
$PLINK --bfile all.clean-inds --geno ${GENO} --maf ${MAF} --hwe ${HWE} --make-bed --out all.clean-base

######################
### WRAP UP ##########
######################

### copy output back to the out path
cp all.clean-base.{bed,bim,fam} ${OUTPATH}
cp missingness.png ancestry* fail-* all.clean-base.log ${OUTPATH}/report/data/
cp ${PBS_O_WORKDIR}/aw_params.sh ${OUTPATH}/report/data/

## generate report
cd ${PBS_O_WORKDIR}
cp -r assets ${OUTPATH}/report/
Rscript -e "library(knitr); knit('report.Rhtml', output='${OUTPATH}/report/report.html')" $PARAMS

### housekeeping
chmod -R 750 $OUTPATH

