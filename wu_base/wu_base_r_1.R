# Julian Q. Zhou
# https://github.com/julianqz

# Install dependent R packages


REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

# for getting versions of installed packages
install.packages("versions")

# Immcantation dependencies in CRAN
# only optparse and BiocManager not in jupyter/datascience-notebook
# c("Rcpp", "devtools", "roxygen2", "testthat", "pkgbuild", "rmarkdown","knitr", "optparse", "BiocManager")
install.packages(c("igraph", "optparse", "BiocManager"))

# Alakazam dependencies in BioConductor
BC_DEP = c("Biostrings", "GenomicAlignments", "IRanges")
BiocManager::install(BC_DEP)

