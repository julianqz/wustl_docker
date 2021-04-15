# Build Docker images for projects at WUSTL

## Cell Ranger (10x Genomics)

`DockerFile_cellranger-6.0.1`

To pull built image: `docker pull julianqz/cellranger-6.0.1:latest`

## Customized image for computational immunology

`DockerFile_JZ_v1`
* Builds on top of `jupyter/datascience-notebook`
* Additional Python libraries: `scanpy`, `presto`, `changeo`
* Additional R packages: `alakazam`, `shazam`, `tigger`, `scoper` (via `install_R_pkgs.R`)
* Versions of the additional libraries/packages are set to be fixed and are reported in `work/docker_build_info.txt`

`install_R_pkgs.R`
* Installs `versions` (for `versions::available:versions`)
* Installs Immcantation dependencies on CRAN that are not part of `jupyter/datascience-notebook`: `optparse`, `BiocManager`
* Installs Immcantation dependencies on BioConductor: `Biostrings`, `GenomicAlignments`, `IRanges`
* During installation of `alakazam` etc., the fixed versions specified in `DockerFile_JZ_v1` are compared against the latest versions on CRAN. If unmatched, archived versions matching the fixed versions are pulled from archive.

Approx. build time on MacBook Pro: 926.7s (896.8s for installation of R packages)

To pull built image: `docker pull julianqz/wustl_jz:v1`
