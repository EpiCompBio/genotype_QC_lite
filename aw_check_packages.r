install = function(x) !x%in%installed.packages()[,"Package"]
if (install("flashpcaR"))
	devtools::install_github("gabraham/flashpca/flashpcaR")
if (install("geneplotter")) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("geneplotter", ask=F)
}
for (p in c("knitr","kableExtra","data.table","gtools","dplyr"))
	if (install(p))
		install.packages(p, quiet=T, dependencies=T, repos='http://cran.us.r-project.org')
