### paths
BASEPATH=/groupvol/med-bio/epiUKB/Airwave/

ZCALLINPUT=true
INFILES=${BASEPATH}/coreExome_zcall/*.zcall
OUTPATH=${BASEPATH}/scripts/genotype_QC_lite/test_output_illumina/
AFFYIDS=
ETHNICFILE=
#ZCALLINPUT=false
#INFILES=${BASEPATH}/affymetrix_genotype/batches_raw/*.bed
#OUTPATH=${BASEPATH}/scripts/genotype_QC_lite/test_output_affy/
#AFFYIDS=${BASEPATH}/affymetrix_genotype/standard/Axiom_UKB_VariantIDmapping.csv
#ETHNICFILE=${BASEPATH}/affymetrix_genotype/ethnicity.txt

STRANDPATH=${BASEPATH}/strandFiles/humancoreexome-12v1-1_a
HIGHLDFILE=${BASEPATH}/../resources/QCfiles/high-LD-regions.txt
HAPMAPPATH=${BASEPATH}/../resources/QCfiles/hapmap3r2_CEU.CHB.JPT.YRI.
HAPMAPBFILE=${HAPMAPPATH}founders.no-at-cg-snps
HAPMAPSNPS=${HAPMAPPATH}no-at-cg-snps.txt

### options
GENOMEBUILD=b37
GENO=0.05
MAF=0.01
HWE=0.00001
HET=3
IMISS=0.03
IBD=0.75
