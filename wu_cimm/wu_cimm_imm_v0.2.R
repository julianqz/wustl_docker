ARGS = commandArgs(trailingOnly=T)

# order of command line arguments should match this
PACKAGES = c("shazam", "tigger", "scoper", "dowser")

stopifnot(length(ARGS)==length(PACKAGES))

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

VERSIONS_WANTED = ARGS
names(VERSIONS_WANTED) = PACKAGES


for (pkg in PACKAGES) {
    
    pkg_ver = VERSIONS_WANTED[pkg]

    #* special treatment 
    #* tigger v1.0.0 (need bug fix for genotypeFasta)
    if (pkg=="tigger" & pkg_ver=="1.0.0") {
        fn = paste0(pkg, "_", pkg_ver, "_fix.tar.gz")
        devtools::install_local(fn)
    } else {
        devtools::install_version(pkg, version=VERSIONS_WANTED[pkg])
    }
    
}
