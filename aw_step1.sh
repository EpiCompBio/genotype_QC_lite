module load plink

PARAMS=/groupvol/med-bio/epiUKB/Airwave/scripts/genotype_QC_lite/aw_params.sh
source $PARAMS
CMDPATH=$(dirname $PARAMS)/

### per-batch preprocessing and QC
cd $CMDPATH
qsub ${CMDPATH}/aw_step2.sh

### create single dataset
# only markers that pass --geno filter in all batches
# only samples that pass missingness and heterogeneity filters
cat *.snplist | LC_ALL=C sort | LC_ALL=C uniq -c | awk -v n="$(ls *.snplist | wc -l)" '$1==n' > n.snplist.all
ls | egrep '^n[0-9]+.bed$' | sed 's/.bed//g' > batches.list
plink --merge-list batches.list --make-bed --out all
plink --bfile all --extract n.snplist.all --make-bed --out all.shared-snps

### create smaller, less redundant dataset
plink --bfile all.shared-snps --indep-pairwise 50 5 0.2 --out all.shared-snps
plink --bfile all.shared-snps --extract all.shared-snps.prune.in --genome --out all.shared-snps

### IBD individuals
Rscript ${PBS_O_WORKDIR}/plot-IBD_modified.R all.shared-snps
# TODO write all.fail-IBD-check.FAILED_QC to fail-IBD-QC.txt?

### divergent ancestry individuals
plink --bfile all.shared-snps --extract ${HAPMAPSNPS} --make-bed --out all.hapmap_snps
plink --bfile all.hapmap_snps --bmerge ${HAPMAPBFILE}.bed ${HAPMAPBFILE}.bim ${HAPMAPBFILE}.fam \
	--extract all.shared-snps.prune.in --make-bed --out all.shared_hapmap_pruned
# TODO Error: 100 variants with 3+ alleles present.
# TODO do PCA

### remove individuals not passing dataset-wide QC


### remove markers not passing dataset-wide QC
plink --bfile XX --maf 0.01 --hwe 0.00001 --out XX
plink --bfile XX --missing --out XX
Rscript lmiss-hist.Rscript
plink --bfile XX --exclude XX --maf 0.01 --hwe 0.00001 --geno 0.05 --make-bed --out XX

