# Julian Q. Zhou
# https://github.com/julianqz

# Install Seurat

ARGS = commandArgs(trailingOnly=T)

# order of command line arguments should match this
PACKAGES = c("anndata", "Seurat")
stopifnot(length(ARGS)==length(PACKAGES))

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

VERSIONS_WANTED = ARGS
names(VERSIONS_WANTED) = PACKAGES


for (pkg in PACKAGES) {

    devtools::install_version(pkg, version=VERSIONS_WANTED[pkg])

}


# for importing .h5ad
#devtools::install_github('satijalab/seurat-data', upgrade="never")
#remotes::install_github("mojaveazure/seurat-disk", upgrade="never")

