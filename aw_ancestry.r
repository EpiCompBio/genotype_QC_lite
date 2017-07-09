library(flashpcaR)
library(MASS)

getContour = function(fx, fy) {
	if (length(fx)>2000) {
		set.seed(1)
		s = sample(length(fx), 2000)
		n = length(s)
		x = fx[s]
		y = fy[s]
	} else {
		n = length(fx)
		x = fx
		y = fy
	}
	den = kde2d(x,y, n=n) 
	z <- array() 
	for (i in 1:n){ 
	        z.x <- max(which(den$x < x[i])) 
	        z.y <- max(which(den$y < y[i])) 
	        z[i] <- den$z[z.x, z.y] 
	} 
	border = quantile(z, probs=0, na.rm=TRUE)
	list(den, border)
}

colMat = t(matrix(c(
"white_British", "#FF280D",
"white_Irish", "#E80C9C",
"white_other", "#9D00FF",
"black_Caribbean", "#0DBEFF",
"black_African", "#00FF1A",
"black_other", "#0CE8A5",
"asian_Indian", "#E8D70C",
"asian_Pakistani", "#FFB500",
"asian_Bangladeshi", "#E8720C",
"Chinese", "#8FFF0D",
"asian_other", "#888888",
"mixed_WhiteAsian", "#888888",
"mixed_WhiteBlackCaribbean", "#888888",
"mixed_WhiteBlackAfrican", "#888888",
"mixed_other", "#888888",
"other", "#888888"
),2))
colnames(colMat) = c("self_ethnic", "col")

mindistThr = as.numeric(as.character(commandArgs(TRUE)[1]))
ethnicityFilename = commandArgs(TRUE)[2]

project_samples = read.table("all.fam",stringsAs=F)[,2]
all_samples_tab = read.table("all.shared_hapmap_pruned.fam",stringsAs=F)[,1:2]
all_samples = all_samples_tab[,2]
wh_project = which(all_samples%in%project_samples)
wh_hapmap = which(!all_samples%in%project_samples)

f  = flashpca("all.shared_hapmap_pruned", ndim=2)
clusters = cutree(hclust(dist(f$vectors[wh_hapmap,])), k=3)
centroids = t(sapply(unique(clusters), function(i) colMeans(f$vectors[wh_hapmap,][clusters==i,])))
x = f$vectors[wh_project,]
distances = t(apply(x, 1, function(i) apply(centroids, 1, function(c) dist(rbind(i,c)))))
wh_cluster = as.numeric(names(tail(sort(table(apply(distances, 1, which.min))),1)))
mindist = distances[,wh_cluster]
summary(mindist)
wh_outlier = wh_project[which(mindist>mindistThr)]
wh_select = wh_project[which(!mindist>mindistThr)]

write.table(all_samples_tab[wh_outlier,], "fail-ancestry-qc.txt", qu=F,ro=F,co=F,sep=" ")

if (!is.na(ethnicityFilename)) {
        png("ancestry.png",1250,800,res=120)
        par(mfrow=c(1,2))
        layout(matrix(c(1,1,2),nrow=1))
	par(mar=c(4,5,1,0))
        plot(f$vectors[,1], f$vectors[,2], col=NULL, xlab="PC1", ylab="PC2", cex.lab=2, cex.axis=1.5, lmain="")
        ethnicityTab = read.table(ethnicityFilename, header=T)
        ethnicityTab = merge(ethnicityTab, colMat, sort=F)
        cols = sapply(project_samples, function(i) as.character(ethnicityTab$col[ethnicityTab$iid==i][1]))
        points(f$vectors[wh_project,1], f$vectors[wh_project,2], col=cols, pch=19)
        points(centroids, col="black", cex=2, pch=19)
        contourData = getContour(f$vectors[wh_select,1], f$vectors[wh_select,2])
        contour(contourData[[1]], levels=contourData[[2]], col="black", add=TRUE, lwd=2)
	par(mar=c(0,0,0,0))
        plot(Inf, axes=F, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
        legend("left", legend=colMat[,1], fill=colMat[,2], cex=1.8, bty="n")
        ethnicity = sapply(project_samples, function(i) as.character(ethnicityTab$self_ethnic[ethnicityTab$iid==i][1]))
        ethnicityTT = table(ethnicity, mindist>mindistThr)
        colnames(ethnicityTT) = c("Kept", "Filtered")
        write.csv(ethnicityTT, "ancestry.csv")
} else {
        png("ancestry.png", 825,800,res=120)
	par(mar=c(4,5,1,0))
        plot(f$vectors[,1], f$vectors[,2], col=NULL, main="", xlab="PC1", ylab="PC2", cex.lab=1.5, cex.axis=1.2)
        points(f$vectors[wh_outlier,1], f$vectors[wh_outlier,2], col="#ff0000ff", pch=19)
        points(f$vectors[wh_select,1], f$vectors[wh_select,2], col="#75c18188", pch=19)
        points(centroids, col="black", pch=19, cex=2)
}
dev.off()

#f10 = flashpca("all.shared_hapmap_pruned", ndim=10)
#f10Mat = f10$vectors[wh_project,]
#rownames(f10Mat) = project_samples
#colnames(f10Mat) = paste0("PC",1:10)
#write.csv(f10Mat, "ancestry_pcs.csv")

