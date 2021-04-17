# https://dowser.readthedocs.io/en/latest/install/

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

# first install dependencies
DEPS = c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp")

for (DEP in DEPS) {
	# if not installed (NA), install
	if (is.na(versions::installed.versions(pkgs=DEP))) {
		cat("Installing dependency:", DEP, "...\n")
		install.packages(DEP)
	}
}

# install the latest development code via devtools, along with devel version of alakazam
library(devtools)
install_bitbucket("kleinstein/alakazam@master")
install_bitbucket("kleinstein/dowser@master")
