args1 = commandArgs(trailingOnly = TRUE)

x = read.table(paste(file.path(getwd(),as.character(args1)), '.lmiss', sep=''), header=T)

png(file.path(getwd(),'missingness.png'), 1200, 1000, res=200)

par(mar=c(4,4.5,0.5,0.5))
hist(log10(x$F_MISS), xlim=c(-4,0), col="#75c181", main="", axes=F, border="#494d55",
        ylab="Number of SNPs", xlab="Fraction of missing data")
abline(v=log10(0.05),lty=2, col="red")
axis(2)
axis(1, at=seq(-4,0), labels=c(0.0001,0.001,0.01,0.1,1))

dev.off()

