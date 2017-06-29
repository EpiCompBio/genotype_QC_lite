args1 = commandArgs(trailingOnly = TRUE)

x = read.table(paste(file.path(getwd(),as.character(args1)), '.lmiss', sep=''), header=T)

png(file.path(getwd(),'missingness.png'))
hist(log10(x$F_MISS), axes=F, xlim=c(-4,0), col="red", ylim=c(0,1e5),
	ylab="Number of SNPs", xlab="Fraction of missing data",main="All SNPs")
axis(side=2,labels=F)
mtext(paste0(seq(0,100,by=20),"K"),side=2,las=2, at=seq(0,1e5,by=2e4),line=1)
axis(side=1,labels=F)
mtext(10**seq(-4,0,by=1),side=1,at=seq(-4,0,by=1),line=1)
abline(v=log10(0.05),lty=2)
dev.off()

