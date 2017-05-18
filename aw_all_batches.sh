#!/bin/sh

#PBS -N all_batches
#PBS -k oe
#PBS -l select=1:ncpus=20:mem=32gb
#PBS -l walltime=8:00:00
#PBS -q med-bio

source ${PBS_O_WORKDIR}/aw_params.sh

module load plink R gcc

######################
### MERGE ############
######################

### create single dataset
# only markers that pass --geno filter in all batches
# only samples that pass missingness and heterogeneity filters
cd ${OUTPATH}
cat *.snplist | LC_ALL=C sort | LC_ALL=C uniq -c | awk -v n="$(ls *.snplist | wc -l)" '$1==n' > ${TMPDIR}/n.snplist.all
ls | egrep '^n[0-9]+.bed$' | sed 's/.bed//g' > ${TMPDIR}/batches.list
plink --merge-list ${TMPDIR}/batches.list --make-bed --out ${TMPDIR}/all
cd ${TMPDIR}
plink --bfile all --extract n.snplist.all --make-bed --out all.shared-snps

### if needed, map Affy IDs to RS IDs
if [[ ! -z "$AFFYIDS" ]];then
  sed -e 's/\"//g' ${AFFYIDS} > mapping.txt
  plink --bfile all.shared-snps --update-name mapping.txt 3 1 20 --make-bed --out all.shared-snps
fi

######################
### QC INDIVIDUALS ###
######################

### create smaller, less redundant dataset
plink --bfile all.shared-snps --exclude ${HIGHLDFILE} --range --indep-pairwise 50 5 0.2 --out all.shared-snps
plink --bfile all.shared-snps --extract all.shared-snps.prune.in --genome --out all.shared-snps

### IBD individuals
Rscript ${PBS_O_WORKDIR}/plot-IBD_modified.R all.shared-snps ${IBD}

### update HapMap data to the same genome build
mkdir -p hapmap
HAPMAPBFILENEW=hapmap/$(basename $HAPMAPBFILE).${GENOMEBUILD}
plink --bfile ${HAPMAPBFILE} --update-alleles ${STRANDPATH}.update_alleles.txt \
      --make-bed --out ${HAPMAPBFILENEW}
${PBS_O_WORKDIR}/update_build.sh ${HAPMAPBFILE} ${STRANDPATH}-${GENOMEBUILD}.strand ${HAPMAPBFILENEW}

### merge with HapMap data
cut -f2 all.shared-snps.bim | sort > all.shared-snps.sort
cut -f2 ${HAPMAPBFILENEW}.bim | sort > ${HAPMAPBFILENEW}.sort
comm -12 all.shared-snps.sort ${HAPMAPBFILENEW}.sort > hapmap_common.snps
plink --bfile all.shared-snps --extract hapmap_common.snps --make-bed --out all.hapmap_snps
plink --bfile all.hapmap_snps --bmerge ${HAPMAPBFILENEW} --out all.missnp
plink --bfile ${HAPMAPBFILENEW} --flip all.missnp.missnp --make-bed --out ${HAPMAPBFILENEW}.flipped
plink --bfile all.hapmap_snps --bmerge ${HAPMAPBFILENEW}.flipped --extract all.shared-snps.prune.in \
      --make-bed --out all.shared_hapmap_pruned

### individuals with divergent ancestry
Rscript ${PBS_O_WORKDIR}/aw_ancestry.r ${ETHNICFILE}

### remove individuals not passing dataset-wide QC
cat fail-* | sort -k1 | uniq | cut -f1,2 -d' '> fail-qc-inds.txt
plink --bfile all.shared-snps --remove fail-qc-inds.txt --make-bed --out all.shared-snps.clean-inds

######################
### QC MARKERS #######
######################

### check missing rate
# markers with missing rates higher than the per-batch threshold are unlikely
plink --bfile all.shared-snps.clean-inds -missing --out all.shared-snps.clean-inds
Rscript ${PBS_O_WORKDIR}/lmiss-hist_modified.R all.shared-snps.clean-inds

### remove markers not passing dataset-wide QC
plink --bfile all.shared-snps.clean-inds --geno ${GENO} --maf ${MAF} --hwe ${HWE} --make-bed --out all.clean-base


### copy output back to the out path
cp all.clean-base.{bed,bim,fam} missingness.png ancestry.* IBD.pdf fail-* n.snplist.all ${OUTPATH}

### compute some filtering stats
Rscript ${PBS_O_WORKDIR}/aw_stats.r ${PBS_O_WORKDIR}/aw_params.sh

