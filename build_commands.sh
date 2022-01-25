#!/usr/bin/env bash


### versions

# base   0.1.0  Ubuntu 20.04.2 LTS (includes g++ 9.3.0)
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

# ref   0.1.x   igblastn 1.17.1
#       0.1.x   IMGT references (standard): release202113-2
#               IMGT references (C57BL6): -
#       0.1.1   IMGT references (standard): release202113-2
#               IMGT references (C57BL6): release202011-3
#       
#       0.2.x   igblastn 1.18.0 (C gene support for Ig seqs)
#               IMGT references (standard): release202150-3
#               IMGT references (C57BL6): release202011-3

# ref_lsf       added LSF env variable for scanpy to work on RIS



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



### wu_*

# Before build, check version in build command is correct
# If wrong and already built, use `docker image tag SOURCE_IMAGE[:TAG] TARGET_IMAGE[:TAG]` to change it

# After build, check /home/jovyan/work/docker_build_info.txt for paths & versions
# e.g.: `docker run --rm -it julianqz/wu_base:main_0.1.0 bash`


cd "/Users/jqz/Dropbox/wustl/code/docker/"

docker build --file wu_base/wu_base_dockerfile --tag julianqz/wu_base:main_0.1.0 ./wu_base

docker push julianqz/wu_base:main_0.1.0

docker build --file wu_presto/wu_presto_dockerfile --tag julianqz/wu_presto:main_0.1.1 \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_presto

docker push julianqz/wu_presto:main_0.1.1

docker build --file wu_cimm/wu_cimm_dockerfile_v0.1 --tag julianqz/wu_cimm:main_0.1.1 \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_cimm

docker push julianqz/wu_cimm:main_0.1.1

# imgt refs

./wu_ref/wu_ref_imgt_get.sh \
	-a "/Users/jqz/Dropbox/common/germline_refs/imgt_zip/" \
	-b "/Users/jqz/Dropbox/common/germline_refs/imgt_select/" \
	-c "/Users/jqz/Dropbox/common/germline_refs/imgt_download.log" \
	-d true

Rscript ./wu_ref/wu_ref_imgt_dedup.R "202113-2" #*

# note context

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_presto:ref_0.1.1 --build-arg BASE_CONTAINER="wu_presto:main_0.1.1" \
	"/Users/jqz/Dropbox/"

docker push julianqz/wu_presto:ref_0.1.1

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_cimm:ref_0.1.1 --build-arg BASE_CONTAINER="wu_cimm:main_0.1.1" \
	"/Users/jqz/Dropbox/"

docker push julianqz/wu_cimm:ref_0.1.1

# add lsf

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_presto:ref_0.1.1_lsf --build-arg BASE_CONTAINER="wu_presto:ref_0.1.1" .

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_cimm:ref_0.1.1_lsf --build-arg BASE_CONTAINER="wu_cimm:ref_0.1.1" .

docker push julianqz/wu_presto:ref_0.1.1_lsf

docker push julianqz/wu_cimm:ref_0.1.1_lsf



### presto 0.1.1d (with debugging for EE set)
# based on wu_base:main_0.1.0 (not 0.1.1 -- 0.1.1 doesn't exist yet)

docker build --file wu_presto/wu_presto_dockerfile_debug --tag julianqz/wu_presto:main_0.1.1d \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_presto

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_presto:ref_0.1.1d --build-arg BASE_CONTAINER="wu_presto:main_0.1.1d" \
	"/Users/jqz/Dropbox/"

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/wu_presto:ref_0.1.1d_lsf --build-arg BASE_CONTAINER="wu_presto:ref_0.1.1d" .

docker push julianqz/wu_presto:ref_0.1.1d_lsf



### ubuntu

cd "/Users/jqz/Dropbox/wustl/code/docker/"

docker build --file DockerFile_ubuntu_main --tag julianqz/ubuntu:main .

docker push julianqz/ubuntu:main

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/ubuntu:main_lsf --build-arg BASE_CONTAINER="ubuntu:main" .

docker push julianqz/ubuntu:main_lsf

# with aspera

cd "/Users/jqz/Dropbox/wustl/code/docker/"

docker build --file DockerFile_ubuntu_aspera --tag julianqz/ubuntu:aspera "/Users/jqz/Dropbox/common/tools"

docker push julianqz/ubuntu:aspera

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/dockerfile_lsf" \
	--tag julianqz/ubuntu:aspera_lsf --build-arg BASE_CONTAINER="ubuntu:aspera" .

docker push julianqz/ubuntu:aspera_lsf

