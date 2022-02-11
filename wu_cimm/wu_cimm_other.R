# Julian Q. Zhou
# https://github.com/julianqz

# Install R packages of choice

PACKAGES = commandArgs(trailingOnly=T)

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

for (pkg in PACKAGES) {
    
    install.packages(pkg)

}
