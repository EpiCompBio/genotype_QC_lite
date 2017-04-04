#############################################
# Outlier analysis using aberrant

# Author: Antonio J Berlanga-Taylor
# Date: 04 April 2017

#Purpose
#=======
# Automatic selection of outliers using 'aberrant' R package
# Here written for genotype analysis but could be run for any other data set with two 
# numerical variables

# For genotype QC, useful to select a homogeneous set of samples to use as reference for marker QC 
# In order to avoid artefacts from population structure. Excluded samples can be re-introduced later.
# PCs from genotype data can be obtained running experiment data against 1000G (or Hapmap) as in UKB appendix 1
# e.g. with MAF >5%, HWE 10^-6, etc for Hapmap or 1000G and then projecting onto these.
# Other values such as heterozygosity, missingness, etc could be used.


#Methods
#=======
# "simple, but robust, algorithm to identify samples with atypical summaries of genome-wide variation"
# See:
# http://www.well.ox.ac.uk/software
# http://www.well.ox.ac.uk/~spencer/Aberrant/aberrant-manual.pdf
# https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/28/1/10.1093/bioinformatics/btr599/2/btr599.pdf?Expires=1491647512&Signature=EglHJtcPkD33-SLCHZYEN6w8L03t8Nqw8OWTlucI-ktjng4HpGHoqyiq-MhuWkghIkmP7yqH8rvThuoRqrS6bGnIGA-rhDQkws2noafo9ERP5TI0meonJ4AubC4olpq1udOc~x5H8Zc55djFuWCkWGtPW~7PhHMGz5T1VPPuZ~9rzRfrBR2aeRzKiGB2ZwTEjLVO8iT8F0UlTtD55XMQ59WaxhzpnLOwmIv0ISXMcgWra6MFiEB4hJ6H11LFqOg4Kn4G7Mf-s7qmRkAr93aVcnIL~XXLaluIMevNirf5NJIP0TTC5p6S7UgxZiL4t8-EFmVt3q8-u~swF4R3W7d45g__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q
# See Fig. 1 for example

# Two numerical variables should be passed.

# Input
# File should be passed without headers and the first two columns with the summary 
# variables.

# TO DO:
# Pass arguments so that user picks col1, col2, col with sample IDs


# Output
# Aberrant plot 
# File with outlier IDs


#Usage
#=====
# To use type:
#    xxx.R [options] [arguments]
#    xxx.R --help


#Options
#=======
#-I    input file name.
#-S    output file name.
#-L    log file name.


# Use docopt, see
#https://github.com/AntonioJBT/various.dir/blob/master/Notes-common-cmds/docopt_argument_parser.txt
#https://github.com/docopt/docopt.R

#############################################


#############################################
# Logging
# TO DO: move to a separate script

##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('~/Desktop/Downloads_to_delete/tests/aberrant_tests/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_aberrant", Sys.Date(), ".txt", sep = ""))
output_file
sink(output_file, append = TRUE, split = TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
# load('/data.dir/R_session_saved_image_aberrant.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_aberrant','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################


#############################
# Import libraries:
library(data.table)

# wget http://www.well.ox.ac.uk/~spencer/Aberrant/aberrant_1.0.tar.gz
# install.packages('aberrant_1.0.tar.gz', repos = NULL, type = "source")
# Gives an error, so untar and create the file 'NAMESPACE' containing the single line
# exportPattern( "." )
# install.packages('aberrant', repos = NULL, type = "source")
library(aberrant)
#############################


#############################
# Set-up arguments:
PC_file <- as.character(args[1])
#############################


#############################
# Example data

# # With flashPCA:
# # https://github.com/gabraham/flashpca/blob/master/flashpcaR/vignettes/flashpcaR.Rmd
# # source("https://bioconductor.org/biocLite.R")
# # devtools::install_github("gabraham/flashpca/flashpcaR")
# library(flashpcaR)
# data(hm3.chr1)
# pca_flash <- scale2(hm3.chr1$bed)
# dim(pca_flash)
# f <- flashpca(X = pca_flash, ndim = 2, stand = 'center')
# PC_file <- as.data.table(f$projection)

# # The following are SNPs, not individuals, but works:
# # http://www.stats.ox.ac.uk/~davison/software/shellfish/snpload-aff.map.gz
# # http://www.stats.ox.ac.uk/~davison/software/shellfish/shellfish.php
# # gunzip snpload-aff.map.gz ; cat snpload-aff.map | cut -f 8,9 > snpload-aff.map.pca
# PC_file <- 'snpload-aff.map.pca'

# # Simulate a dataset, the flashPCA Hapmap3 has no outliers
# N <- 10000
# M <- 2
# sim1 <- matrix(rnorm(N * M, mean = 0, sd = 1), N, M)
# sim1
# PC_file <- as.data.table(sim1)

# # Simulate a dataset with two clusters:
# N <- 500
# M <- 2
# sim2 <- matrix(rnorm(N * M, mean = 4, sd = 1), N, M)
# sim_2clust <- rbind(sim1, sim2)
# plot(sim_2clust)
# PC_file <- as.data.table(sim_2clust)
#############################

# Read files:
PC_file <- fread(PC_file, header = F, stringsAsFactors = FALSE, select = c(1:2))
PC_file
head(PC_file)
dim(PC_file)
tail(PC_file)
summary(PC_file)
#############################


#############################
# Run aberrant
# Use summary statistics, and/or: missingness, ancestry, probe intensity, gender separately
# aberrant with lambda set to 20 for ancestry PC1 and PC2 as summary stats:
aberrant_1 <- aberrant(x = PC_file, lambda = 20)
class(aberrant_1)
str(aberrant_1)

# Write outliers to file:
fwrite(as.data.table(aberrant_1$outlier), 'aberrant_outlier_IDs.txt', 
       sep = '\t', col.names = F)

# Plot:
pdf('aberrant_plot.pdf', width = 8, height = 5)
aberrant.plot(aberrant_1)
dev.off()

# Legend from original paper:
# Fig. 1.
"Outlier identification for xxx samples genotyped on xxx. 
‘Normal’ individuals are coloured from black to grey, darker colours denoting higher 
density of individuals. 
Outliers are coloured from orange to red, with redder colours denoting 
higher posterior probability of being an outlier. 
The 99% confidence ellipse of the inferred distribution of ‘normal’ individuals is 
shown as a dashed grey line."
#############################


#############################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx
#############################