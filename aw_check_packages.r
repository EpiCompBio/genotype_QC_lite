install = function(x) !x%in%installed.packages()[,"Package"]
if (install("flashpcaR"))
	devtools::install_github("gabraham/flashpca/flashpcaR")
if (install("geneplotter")) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("geneplotter", ask=F)
}

