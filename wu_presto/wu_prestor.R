# Julian Q. Zhou
# https://github.com/julianqz

# Install prestor

# https://bitbucket.org/kleinstein/prestor/src/master/README.md

REPO_CRAN = "http://cran.us.r-project.org"
#REPO_CRAN = "http://lib.stat.cmu.edu/R/CRAN/"

# need this to install from within Docker
options(repos=REPO_CRAN)

# first install dependencies
#* already in jupyter/datascience-notebook
#install.packages(c("devtools", "roxygen2"))

# install the latest development code via devtools
#* "javh/prototype-prestor@default" doesn't work any more -- repo has changed location
library(devtools)
install_bitbucket("kleinstein/prestor@master")
