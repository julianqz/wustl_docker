# Docker containers for bioinformatics work at WUSTL

## Ubuntu v20.04

* Minimal Ubuntu with essentials like `git`, `vim`, `lftp` etc. 

`docker pull julianqz/ubuntu:main` (102 MB) 

* Same as above, plus Aspera Connect v4.1.1.73 

`docker pull julianqz/ubuntu:aspera` (322 MB)

## Cell Ranger v6.0.1

Note: Reference genome assemblies are not included and need to be downloaded separately.

`docker pull julianqz/cellranger-6.0.1` (1.1 GB)

## Customized containers for computational immunology

See [`build_commands.sh`](./build_commands.sh) for version info. 

* v0.1.1 is the main container version being used. 
* v0.2.0 will be used for future projects.

Within each container, installed versions are recorded at `work/docker_build_info.txt`.

For general use by the public, use the pull commands listed under [`4)`](https://github.com/julianqz/wustl_docker#4-presto--cimm-container-with-igblast-and-igblast-ready-imgt-ig-references).

### 1) Base container

* [`jupyter/datascience-notebook`](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html): Ubuntu, Python, R
* R packages required for dependency: alakazam, knitr, BiocManager, etc.

### 2) pRESTO container (for preprocessing bulk BCR sequencing)

* Base container
* pRESTO

This container is designed to work with pipeline scripts in [`wu_presto/pipelines`](./wu_presto/pipelines) and helper scripts in [`wu_presto/scripts`](./wu_presto/scripts), some of which were adapted from [Immcantation](https://bitbucket.org/kleinstein/immcantation). 

Together, they support preprocessing sequencing libraries prepared using the NEB NextImmune kits. 

Primer sequences for the human kit were from [Immcantation](https://bitbucket.org/kleinstein/immcantation) ([`wu_presto/protocols/AbSeq`](./wu_presto/protocols/AbSeq)). 

Primer sequences for the mouse kit were from NEB's [Galaxy workflow](https://usegalaxy.org/u/bradlanghorst/w/presto-nebnext-immune-seq-workflow-v320) (accessed 2021-11-18) with renaming ([`wu_presto/protocols/NEB_NextImmune_mouse`](./wu_presto/protocols/NEB_NextImmune_mouse)).

### 3) `cimm` (<u>c</u>omputational <u>i</u>mmunology) container (BCR & GEX processing & analysis)

* Base container
* BCR parsing: Change-O
* BCR analysis: TIgGER, SHazaM, SCOPer
* Phylogenetic inference: PHYLIP, IgPhyML, (dowser)
* Single-cell gene expression: Scanpy, Seurat, (scRepertoire)
* GenBank submission: (tbl2asn)

Packages in parentheses are not installed in all versions of the container.

### 4) pRESTO / `cimm` container with IgBLAST and IgBLAST-ready IMGT Ig references

* pRESTO / `cimm` container
* IgBLAST
* IMGT Ig references
* IgBLAST database files built based on IMGT Ig references

`docker pull julianqz/wu_presto:ref_0.1.1` (1.7 GB)

`docker pull julianqz/wu_presto:ref_0.2.0` (1.8 GB)

`docker pull julianqz/wu_cimm:ref_0.1.1` (2.2 GB)

`docker pull julianqz/wu_cimm:ref_0.2.0` (2.1 GB)

### 5) Containers tailored to run on the LSF system on WUSTL RIS

* pRESTO / `cimm` container with IgBLAST and IgBLAST-ready IMGT Ig references
* Added environment variable `NUMBA_CACHE_DIR`
* Added `ENTRYPOINT`

`docker pull julianqz/wu_presto:ref_0.1.1_lsf`

`docker pull julianqz/wu_presto:ref_0.2.0_lsf`

`docker pull julianqz/wu_cimm:ref_0.1.1_lsf`

`docker pull julianqz/wu_cimm:ref_0.2.0_lsf`

## Customized container for reproducing published analysis

* R v4.1.0, igraph v1.2.5, ggplot2 v3.3.5, ggtree v3.0.4, circlize v0.4.13

`docker pull julianqz/wu_pub:r_4.1.0` (1.0 GB)

