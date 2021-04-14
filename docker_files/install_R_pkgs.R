ARGS = commandArgs(trailingOnly=T)
stopifnot(length(ARGS)==4)

# order of command line arguments should match this
PACKAGES = c("alakazam", "shazam", "tigger", "scoper")

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

VERSIONS_WANTED = ARGS
names(VERSIONS_WANTED) = PACKAGES

# latest version on CRAN
install.packages("versions")

VERSIONS_LATEST = sapply(PACKAGES, function(s){
    return(versions::available.versions(s)[[1]][["version"]][1])
}, USE.NAMES=T)

# Immcantation dependencies in CRAN
# only optparse and BiocManager not in jupyter/datascience-notebook
# c("Rcpp", "devtools", "roxygen2", "testthat", "pkgbuild", "rmarkdown","knitr", "optparse", "BiocManager")
install.packages(c("optparse", "BiocManager"))

# Alakazam dependencies in BioConductor
BC_DEP = c("Biostrings", "GenomicAlignments", "IRanges")
BiocManager::install(BC_DEP)

for (pkg in PACKAGES) {
    if (VERSIONS_LATEST[pkg] == VERSIONS_WANTED[pkg]) {
        # if latest on CRAN is the same as the specified version
        # install latest from CRAN
        install.packages(pkg, clean=T)
    } else {
        # else, install the specified version from archive
        devtools::install_version(pkg,
                                  version=VERSIONS_WANTED[pkg], 
                                  repos=REPO_CRAN)
    }
}
