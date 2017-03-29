module load plink
module load R
# TODO check it's installed: devtools::install_github("gabraham/flashpca/flashpcaR")

set -eu

PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/genotype_QC_lite/aw_params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/

######################
### BATCH QC & MERGE #
######################

### per-batch preprocessing and QC
cd ${CMDPATH}
qsub ${CMDPATH}/aw_step2.sh

cd ${PLINKPATH}

### create single dataset
# only markers that pass --geno filter in all batches
# only samples that pass missingness and heterogeneity filters
cat *.snplist | LC_ALL=C sort | LC_ALL=C uniq -c | awk -v n="$(ls *.snplist | wc -l)" '$1==n' > n.snplist.all
ls | egrep '^n[0-9]+.bed$' | sed 's/.bed//g' > batches.list
plink --merge-list batches.list --make-bed --out all
plink --bfile all --extract n.snplist.all --make-bed --out all.shared-snps

######################
### QC INDIVIDUALS ###
######################

### create smaller, less redundant dataset
plink --bfile all.shared-snps --exclude ${HIGHLDFILE} --range --indep-pairwise 50 5 0.2 --out all.shared-snps
plink --bfile all.shared-snps --extract all.shared-snps.prune.in --genome --out all.shared-snps

### IBD individuals
Rscript ${PBS_O_WORKDIR}/plot-IBD_modified.R all.shared-snps
# TODO write all.fail-IBD-check.FAILED_QC to fail-IBD-QC.txt?

### update HapMap data to the same genome build
mkdir -p hapmap
HAPMAPBFILENEW=${PLINKPATH}/hapmap/$(basename $HAPMAPBFILE).${GENOMEBUILD}
plink --bfile ${HAPMAPBFILE} --update-alleles ${STRANDPATH}.update_alleles.txt \
      --make-bed --out ${HAPMAPBFILENEW}
${PBS_O_WORKDIR}/update_build.sh ${HAPMAPBFILE} ${STRANDPATH}-${GENOMEBUILD}.strand ${HAPMAPBFILENEW}

### merge with HapMap data
plink --bfile all.shared-snps --extract ${HAPMAPSNPS} --make-bed --out all.hapmap_snps
plink --bfile all.hapmap_snps --bmerge ${HAPMAPBFILENEW} --out all.missnp
plink --bfile ${HAPMAPBFILENEW} --flip all.missnp.missnp --make-bed --out ${HAPMAPBFILENEW}.flipped
plink --bfile all.hapmap_snps --bmerge ${HAPMAPBFILENEW}.flipped --extract all.shared-snps.prune.in \
      --make-bed --out all.shared_hapmap_pruned

### individuals with divergent ancestry
Rscript ${PBS_O_WORKDIR}/aw_ancestry.r

### remove individuals not passing dataset-wide QC
cat fail-* | sort -k1 | uniq | cut -f1,2 -d' '> fail-qc-inds.txt
plink --bfile all.shared-snps --remove fail-qc-inds.txt --make-bed --out all.shared-snps.clean-inds

######################
### QC MARKERS #######
######################

### check missing rate
# markers with missing rates higher than the per-batch threshold are unlikely
plink --bfile all.shared-snps.clean-inds -missing --out all.shared-snps.clean-inds
Rscript ${CMDPATH}/lmiss-hist_modified.R all.shared-snps.clean-inds

### remove markers not passing dataset-wide QC
plink --bfile all.shared-snps.clean-inds --geno 0.05 --maf 0.01 --hwe 0.00001 --make-bed --out all.clean-base

