library(flashpcaR)

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
wh_outlier = wh_project[which(mindist>0.005)]

write.table(all_samples_tab[wh_outlier,], "fail-ancestry-qc.txt", qu=F,ro=F,co=F,sep=" ")

png("ancestry.png")
plot(f$vectors[,1], f$vectors[,2], col=NULL, main="PCA with HapMap", xlab="PC1", ylab="PC2")
points(f$vectors[wh_project,1], f$vectors[wh_project,2], col="#00ff0088", pch=19)
points(f$vectors[wh_hapmap,1], f$vectors[wh_hapmap,2], col="#00000008", pch=19)
#points(centroids, col="green", pch=19)
points(f$vectors[wh_outlier,1], f$vectors[wh_outlier,2], col="#ff0000ff", pch=19)
dev.off()

