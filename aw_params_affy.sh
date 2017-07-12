PROJECTNAME="Airwave (Affymetrix)"

### paths
BASEPATH=/groupvol/med-bio/epiUKB/Airwave/

ZCALLINPUT=false
INFILES=${BASEPATH}/affymetrix_genotype/ICL_raw_1-9.bed
OUTPATH=${BASEPATH}/QC/affy_b37_maf0.01/
AFFYIDS=${BASEPATH}/affymetrix_genotype/standard/Axiom_UKB_VariantIDmapping.csv
ETHNICFILE=${BASEPATH}/affymetrix_genotype/ethnicity.txt
SEXFILE=${BASEPATH}/affymetrix_genotype/sex.txt

STRANDPATH=${BASEPATH}/strandFiles/humancoreexome-12v1-1_a
HIGHLDFILE=${BASEPATH}/../resources/QCfiles/high-LD-regions.txt
HAPMAPPATH=${BASEPATH}/../resources/QCfiles/hapmap3r2_CEU.CHB.JPT.YRI.
HAPMAPBFILE=${HAPMAPPATH}founders.no-at-cg-snps
HAPMAPSNPS=${HAPMAPPATH}no-at-cg-snps.txt

### options
GENOMEBUILD=b37
GENO=0.02
MAF=0.01
HWE=0.00001
HET=3
IMISS=0.03
REL=0.1875
ANC=0.005
