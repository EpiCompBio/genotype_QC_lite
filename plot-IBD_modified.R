#################################
#Identity by descent plot
#Plots IBD results from genotype data

#Script modified from www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf

#Requires:
  #Input file with data to read in (.genome file) after plink analysis
  #Set cut-off of IBD, 'IBD_cutoff', currently at 0.1875, typical value used though
  #Runs in current directory

#Antonio Berlanga-Taylor
#David Mosen-Ansorena
#April 2017

#################################

library("data.table")

args = commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = as.character(args[1])
input_path = file.path(getwd(), input_file)
output_path = file.path(getwd(), 'IBD.pdf')
IBD_cutoff = as.numeric(args[2])
FAILED_IDs_file = file.path(getwd(), 'fail-IBD-qc.txt')

#Read in data:
data=fread(paste(input_path, ".genome", sep=''))

#Plot:
pdf(output_path)

hist(data$PI_HAT, ylim=c(0,nrow(data)), col="RED", breaks=100, 
     xlab="Estimated mean pairwise IBD",main="")

abline(v=IBD_cutoff, col="gray32", lty=2)

dev.off()

#Write to table individuals failing IBD results:
out = which(data$PI_HAT > IBD_cutoff)
write.table(data[out,], file = FAILED_IDs_file, append = FALSE, col.names=TRUE, row.names=FALSE, 
            sep='\t', quote= FALSE)


#pihats = seq(0.1,0.8,0.05)
#nfilt = sapply(pihats, function(i) length(unique(data[which(data$PI_HAT > i),FID1])))
#plot(pihats, nfilt, xlab="IBD", ylab="# samples")
