# Bioconductor release version specified via command line
ARGS = commandArgs(trailingOnly=T)
stopifnot(length(ARGS)==1)

# Bioconductor pkg to be installed
PACKAGES = c("scRepertoire")

for (pkg in PACKAGES) {

    # `update` and `ask` deal with prompt to update installed packges
    
    # `checkBuilt` if TRUE considers a pkg built under an earlier major.minor version of R to be old
    
    # `force` if TRUE re-downloads a pkg that is currently up-to-date

    # `version` specifies the Bioconductor release version from which to install pkg 
    #     Defaults to BiocManager::version() ==> 
    #     No need to specify as long as it's made sure that the desired pkg version is the same as 
    #        the latest pkg version available in the Bioconductor release version specified

    BiocManager::install(pkg, version=ARGS,
                         update=FALSE, ask=FALSE, 
                         checkBuilt=TRUE, force=FALSE)

}
