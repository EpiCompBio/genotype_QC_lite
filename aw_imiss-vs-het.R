#Original script: www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf

library("geneplotter")

args = commandArgs(trailingOnly=T)

#Read data:
input_path = file.path(getwd(), as.character(args[1]))
imiss = read.table(paste(input_path, ".imiss", sep=""), h=T, stringsAs=F)
het = read.table(paste(input_path, ".het", sep=""), h=T, stringsAs=F)
het$meanHet = (het$N.NM. - het$O.HOM.) / het$N.NM.

#Set filter thresholds:
HET_SDS = as.numeric(as.character(args[2]))
HET_SD_THRS = mean(het$meanHet,na.rm=T) + c(-HET_SDS,HET_SDS) * sd(het$meanHet,na.rm=T)
IMISS_THR = as.numeric(as.character(args[3]))

#Plot:
colors = densCols(log10(imiss[,6]), het$meanHet)
pdf(file.path(getwd(), 'imiss-vs-het.pdf'))
plot(log10(imiss[,6]), het$meanHet, col=colors, xlim=c(-3,0),
     ylim=c(0,0.5), pch=20, xlab="Proportion of missing genotypes",
     ylab="Heterozygosity rate", axes=F)
axis(2, at=seq(0,0.5,by=0.05), tick=T)
axis(1, at=-3:0, labels=10**c(-3:0))
abline(h=HET_SD_THRS, col="red", lty=2)
abline(v=log10(IMISS_THR), col="RED", lty=2)
dev.off()

#Write failed samples:
prepare = function(df, fail) {
  names(df)[3] = "VALUE"
  if (nrow(df)) {
    cbind(df,FAILED=fail)
  } else {
    df$FAILED=character(0)
    df
  }
}
high_het = het[which(het$meanHet>HET_SD_THRS[2]), c("FID","IID","meanHet")]
low_het = het[which(het$meanHet<HET_SD_THRS[1]), c("FID","IID","meanHet")]
high_imiss = imiss[which(imiss[,6]>IMISS_THR), c("FID","IID","F_MISS")]
failed = rbind(prepare(high_het,"high_het"),prepare(low_het,"low_het"),prepare(high_imiss,"high_imiss"))
write.table(failed, file = file.path(getwd(),'fail-imisshet-qc.txt'), co=T, qu=F, ro=F, sep='\t')

