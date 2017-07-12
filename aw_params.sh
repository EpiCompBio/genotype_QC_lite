PROJECTNAME="Airwave (Illumina)"

### paths
BASEPATH=/groupvol/med-bio/epiUKB/Airwave/

ZCALLINPUT=true
# same as ${BASEPATH}/coreExome_zcall/*.zcall but properly ordered
INFILES=$(ls ${BASEPATH}/coreExome_zcall/ | egrep '^n[0-9]+.zcall$' | sort -n -k1.2 | \
        sed "s@^@${BASEPATH}\/coreExome_zcall\/@" | sed ':a;N;$!ba;s/\n/ /g')
OUTPATH=${BASEPATH}/QC/illumina_b37_maf0.01/
AFFYIDS=
ETHNICFILE=${BASEPATH}/coreExome_genotype/ethnicity.txt
SEXFILE=${BASEPATH}/coreExome_genotype/sex.txt

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
