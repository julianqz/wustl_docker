# Julian Q. Zhou
# https://github.com/julianqz

# Install captioner (removed from CRAN 2023-08-18)

#ARGS = commandArgs(trailingOnly=T)
#stopifnot(length(ARGS)==1)

# order of command line arguments should match this
PACKAGES = c("captioner")
PATHS = c("https://cran.r-project.org/src/contrib/Archive/captioner/captioner_2.2.3.tar.gz")
names(PATHS) = PACKAGES

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

#VERSIONS_WANTED = ARGS
#names(VERSIONS_WANTED) = PACKAGES


for (pkg in PACKAGES) {

    #devtools::install_version(pkg, version=VERSIONS_WANTED[pkg])
    install.packages(PATHS[pkg], type = "source")

}
