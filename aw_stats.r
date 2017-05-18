getLines = function(thepath, thefiles) {
	cmd = paste0("wc -l `ls -v ",thepath,"/",thefiles,"`| head -n -1")
	out = sapply(trimws(system(cmd,intern=T)), function(i) strsplit(i," ")[[1]][1])
	as.numeric(as.character(as.vector(out)))
}
getLines1 = function(thepath, thefile)
	as.numeric(system(paste0("cat ",thepath,"/",thefile," | wc -l"), intern=T))

params = commandArgs(trailingOnly=TRUE)[1]
infiles = system(paste0("source ",params,"; echo $INFILES"), intern=T)
infiles = gsub(".bed","", basename(strsplit(infiles," ")[[1]]))
thepath = system(paste0("source ",params,"; echo $OUTPATH"), intern=T)

df = data.frame(nIndivKept=getLines(thepath,"n*.fam"),
				nIndivFilt=getLines(thepath,"n*.fail-imisshet-qc.txt"),
				nMarkTotal=getLines(thepath,"n*.bim"),
				nMarkKept=getLines(thepath,"n*.snplist"))
df$nIndivTotal = df$nIndivKept + df$nIndivFilt
df$nMarkFilt = df$nMarkTotal - df$nMarkKept
df = df[,c("nIndivTotal","nIndivKept","nIndivFilt","nMarkTotal","nMarkKept","nMarkFilt")]
rownames(df) = infiles

nMarkKept = getLines1(thepath,"n.snplist.all")
totalBatchFilt = c(apply(df[,grep("Indiv",colnames(df))], 2, sum),
					 df[1,"nMarkTotal"], nMarkKept, df[1,"nMarkTotal"]-nMarkKept)
df = rbind(df, totalBatchFilt=totalBatchFilt)
	
jointIndivKept = getLines1(thepath,"all.clean-base.fam")
jointMarkKept = getLines1(thepath,"all.clean-base.bim")
afterJointFilt = c(df["totalBatchFilt","nIndivTotal"], jointIndivKept, 
				   df["totalBatchFilt","nIndivTotal"]-jointIndivKept,
				   df["totalBatchFilt","nMarkTotal"], jointMarkKept, 
				   df["totalBatchFilt","nMarkTotal"]-jointMarkKept)
df = rbind(df, afterJointFilt=afterJointFilt)

write.csv(df, file.path(thepath,"filter_stats.csv"))

