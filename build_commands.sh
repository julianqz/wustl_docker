#!/usr/bin/env bash

# Julian Q. Zhou
# https://github.com/julianqz

# Documentation of build commands used to build Docker containers


### versions

# base    0.1.0 Ubuntu 20.04.2 LTS (includes g++ 9.3.0)
#               Python 3.8.8
#               R 4.0.3
#               alakazam 1.1.0
#               knitr 1.31
#               BiocManager 1.30.12 (Bioconductor 3.12)
#
#         0.2.0 Ubuntu 20.04.3 LTS (no g++)
#               Python 3.9.7
#               R 4.1.1
#               alakazam 1.2.0
#               knitr 1.32
#               BiocManager 1.30.16 (Bioconductor 3.14)
#
#         0.3.0 Ubuntu 22.04.2 LTS [to check, `cat /etc/issue`]
#               Python 3.11.4
#               R 4.3.1
#               alakazam 1.3.0
#               knitr 1.45
#               BiocManager 1.30.22 (Bioconductor 3.18)
#               igraph 2.0.2

# presto 0.1.0  initial
#               pRESTO 0.6.2
#        0.1.1  increase cd-hit-est memory limit
#        0.1.1d debugging for EE set (no fix) 
#
#        0.2.0  pRESTO 0.7.0 (includes parameter to increase cd-hit-est memory limit)
#               prestor fixed compatbility issue with knitr 1.32

# cimm   0.1.0  initial
#               changeo 1.0.2
#               scanpy 1.7.2
#               Seurat 4.0.1
#               PHYLIP 3.697
#               IgPHYML 1.1.3
#               shazam 1.0.2
#               tigger 1.0.0
#               scoper 1.1.0
#               tbl2asn 25.8
#               
#        0.1.1  bug fixes for shazam::distToNearest (v1.0.2) and 
#                             tigger::genotypeFasta (v1.0.0)
#               scikit-misc 0.1.4 (new addition)
#
#
#        0.2.0  changeo 1.2.0
#               scanpy 1.8.2
#               Seurat 4.0.6
#               scRepertoire 1.4.0 via Bioconductor 3.14 (new addition)
#               shazam 1.1.0
#               scoper 1.2.0
#               dowser 0.1.0 via CRAN (new addition)
#               tbl2asn - (removed)
#
#       0.2.2   SeuratData 0.2.2
#               SeuratDisk 0.0.0.9020
#
#       0.3.0   presto 0.7.2
#               changeo 1.3.0
#               scanpy 1.9.8
#               anndata (python) 0.10.5.post1
#               scikit-misc 0.3.1
#               SeuratObject 5.0.1
#               Seurat 5.0.1
#               scRepertoire (removed)
#               tigger 1.1.0
#               shazam 1.2.0
#               scoper 1.3.0
#               dowswer 2.1.0
#               anndata (R) 0.7.5.6
#               circlize 0.4.16
#               IgPHYML 2.0.0
#       0.3.1   muon 0.1.6
#       0.3.2   Matrix.utils 0.9.8
#               pheatmap 1.0.12
#               SingleCellExperiment 1.24.0
#               edgeR 4.0.16
#               DESeq2 1.42.1
#               apeglm 1.24.0

# ref   0.1.x   igblastn 1.17.1
#       0.1.0   IMGT references (standard): release202113-2
#               IMGT references (C57BL6): -
#       0.1.1   IMGT references (standard): release202113-2
#               IMGT references (C57BL6): release202011-3
#       
#       0.2.x   igblastn 1.18.0 (C gene support for Ig seqs)
#               IMGT references (standard): release202150-3
#               IMGT references (C57BL6): release202011-3
#  
#       0.3.x  igblastn 1.22.0
#              IMGT references (standard): release202405-2
#              IMGT references (C57BL6): release202011-3

# ref_lsf       added LSF env variable for scanpy to work on RIS

# beta 0.2.0    (based on wu_base:main_0.2.0)
#               R - anndata, 0.7.5.6
#               R - Seurat, 5.0.0
#               Python - anndata, 0.10.4

# pub  r_4.1.0  R 4.1.0
#               BioConductor 3.13
#               igraph 1.2.5
#               circlize 0.4.13
#               ggtree 3.0.4
#               vioplot 0.3.6
#               ggplot2 3.3.5
#               alakazam 1.1.0
#               shazam 1.0.2 with bug fix for shazam::distToNearest
#               tigger 1.0.0 with bug fix for tigger::genotypeFasta


### design

# jupyter/datascience-notebook
#	
# 	 wu_base:main_: Bioconductor dependencies, alakazam, knitr
#
#		     wu_presto:main_: presto, prestor (alakazam-dependent), muscle, vsearch, CD-HIT, BLAST
#	                          abseq pipeline
#
#	               wu_presto:ref_: IgBLAST, IMGT
#
#	        			     _lsf: ENTRYPOINT
#
#	         wu_cimm:main_: changeo, shazam, tigger, scoper
#		                    scanpy, Seurat, (scRepertoire)
#		                    PHYLIP, igphyml, (dowser)
#	        	            (tbl2asn)
#
#		            wu_cimm:ref_: IgBLAST, IMGT
#
#		    		 	    _lsf: ENTRYPOINT


# jupyter/r-notebook
#
#    wu_pub:r_4.1.0: specific versions of R packages
#           _lsf: ENTRYPOINT


### system clean up: https://docs.docker.com/reference/cli/docker/system/prune/

### wu_base/presto/cimm

# Before build, check version in build command is correct
# If wrong and already built, use `docker image tag SOURCE_IMAGE[:TAG] TARGET_IMAGE[:TAG]` to change it

# After build, check /home/jovyan/work/docker_build_info.txt for paths & versions
# e.g.: `docker run --rm -it julianqz/wu_base:main_0.1.0 bash`


# use --progress=plain with `docker build` for non-truncated output

cd "/Users/jqz/Dropbox/wustl/code/docker/"

docker build --progress=plain --file wu_base/wu_base_dockerfile --tag julianqz/wu_base:main_0.1.0 ./wu_base

docker push julianqz/wu_base:main_0.1.0

docker build --progress=plain --file wu_presto/wu_presto_dockerfile --tag julianqz/wu_presto:main_0.1.1 \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_presto

docker push julianqz/wu_presto:main_0.1.1

docker build --progress=plain --file wu_cimm/wu_cimm_dockerfile_v0.3 --tag julianqz/wu_cimm:main_0.3.2 \
	--build-arg BASE_CONTAINER="wu_base:main_0.3.0" ./wu_cimm

docker push julianqz/wu_cimm:main_0.3.2

# imgt refs

./wu_ref/wu_ref_imgt_get.sh \
	-a "/Users/jqz/Dropbox/common/germline_refs/imgt_zip/" \
	-b "/Users/jqz/Dropbox/common/germline_refs/imgt_select/" \
	-c "/Users/jqz/Dropbox/common/germline_refs/imgt_download.log" \
	-d true

Rscript ./wu_ref/wu_ref_imgt_dedup.R "202113-2" #*

# note context

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_presto:ref_0.1.1 --build-arg BASE_CONTAINER="wu_presto:main_0.1.1" \
	"/Users/jqz/Dropbox/"

docker push julianqz/wu_presto:ref_0.1.1

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_cimm:ref_0.3.2 --build-arg BASE_CONTAINER="wu_cimm:main_0.3.2" \
	"/Users/jqz/Dropbox/"

docker push julianqz/wu_cimm:ref_0.3.2

# add lsf

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_presto:ref_0.1.1_lsf --build-arg BASE_CONTAINER="wu_presto:ref_0.1.1" .

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_cimm:ref_0.3.2_lsf --build-arg BASE_CONTAINER="wu_cimm:ref_0.3.2" .

docker push julianqz/wu_presto:ref_0.1.1_lsf

docker push julianqz/wu_cimm:ref_0.3.2_lsf



### presto 0.1.1d (with debugging for EE set)
# based on wu_base:main_0.1.0 (not 0.1.1 -- 0.1.1 doesn't exist yet)

docker build --progress=plain --file wu_presto/wu_presto_dockerfile_debug --tag julianqz/wu_presto:main_0.1.1d \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_presto

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_presto:ref_0.1.1d --build-arg BASE_CONTAINER="wu_presto:main_0.1.1d" \
	"/Users/jqz/Dropbox/"

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_presto:ref_0.1.1d_lsf --build-arg BASE_CONTAINER="wu_presto:ref_0.1.1d" .

docker push julianqz/wu_presto:ref_0.1.1d_lsf



### ubuntu

cd "/Users/jqz/Dropbox/wustl/code/docker/"

docker build --progress=plain --file dockerfile_ubuntu_main --tag julianqz/ubuntu:main .

docker push julianqz/ubuntu:main

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/ubuntu:main_lsf --build-arg BASE_CONTAINER="ubuntu:main" .

docker push julianqz/ubuntu:main_lsf

# with aspera

cd "/Users/jqz/Dropbox/wustl/code/docker/"

docker build --progress=plain --file dockerfile_ubuntu_aspera --tag julianqz/ubuntu:aspera "/Users/jqz/Dropbox/common/tools"

docker push julianqz/ubuntu:aspera

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/ubuntu:aspera_lsf --build-arg BASE_CONTAINER="ubuntu:aspera" .

docker push julianqz/ubuntu:aspera_lsf



### wu_pub

# use ./wu_cimm for context (location of tar.gz containing bug fixes for shazam & tigger)

docker build --progress=plain --file wu_pub/wu_pub_r_4.1.0 --tag julianqz/wu_pub:r_4.1.0 ./wu_cimm

docker push julianqz/wu_pub:r_4.1.0

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_pub:r_4.1.0_lsf --build-arg BASE_CONTAINER="wu_pub:r_4.1.0" .

docker push julianqz/wu_pub:r_4.1.0_lsf


### beta

docker build --progress=plain --file wu_beta/wu_beta_dockerfile_v0.1 --tag julianqz/wu_beta:main_0.2.0 \
	--build-arg BASE_CONTAINER="wu_base:main_0.2.0" ./wu_beta

docker push julianqz/wu_beta:main_0.2.0

docker build --progress=plain --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_beta:main_0.2.0_lsf --build-arg BASE_CONTAINER="wu_beta:main_0.2.0" .

docker push julianqz/wu_beta:main_0.2.0_lsf

