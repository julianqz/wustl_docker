ARGS = commandArgs(trailingOnly=T)
stopifnot(length(ARGS)==3)

# order of command line arguments should match this
PACKAGES = c("shazam", "tigger", "scoper")

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

VERSIONS_WANTED = ARGS
names(VERSIONS_WANTED) = PACKAGES


for (pkg in PACKAGES) {
    
    #* special treatment 
    #* shazam v1.0.2 (need bug fix for distToNearest)
    #* tigger v1.0.0 (need bug fix for genotypeFasta)
    if (pkg=="shazam" | pkg=="tigger") {
        fn = paste0(pkg, "_", VERSIONS_WANTED[pkg], "_fix.tar.gz")
        devtools::install_local(fn)
    } else {
        devtools::install_version(pkg, version=VERSIONS_WANTED[pkg])
    }
    
}
