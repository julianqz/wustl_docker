# Julian Q. Zhou
# https://github.com/julianqz

# Install knitr

ARGS = commandArgs(trailingOnly=T)
stopifnot(length(ARGS)==1)

# order of command line arguments should match this
PACKAGES = c("knitr")

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

VERSIONS_WANTED = ARGS
names(VERSIONS_WANTED) = PACKAGES


for (pkg in PACKAGES) {

    devtools::install_version(pkg, version=VERSIONS_WANTED[pkg])

}
