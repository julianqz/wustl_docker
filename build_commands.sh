#!/usr/bin/env bash


### versions

# presto 0.1.0  initial
#        0.1.1  increase cd-hit-est limit
#        0.1.1d debugging for EE set (no fix) 

# cimm   0.1.0  initial
#        0.1.1  bug fixes for shazam::distToNearest (v1.0.2) and 
#                             tigger::genotypeFasta (v1.0.0)



### design

# jupyter/datascience-notebook
#	
# 	 wu_base:main_: bioconductor dependencies, alakazam
#
#		     wu_presto:main_: presto, prestor (alakazam-dependent), muscle, vsearch, CD-HIT, BLAST
#	                          abseq pipeline
#
#	               wu_presto:ref_: IgBLAST, IMGT
#
#	        			     _lsf: ENTRYPOINT
#
#	         wu_cimm:main_: changeo, shazam, tigger, scoper
#		                    scanpy, seurat
#		                    PHYLIP, igphyml, (dowser)
#	        	            tbl2asn
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

docker build --file wu_presto/wu_presto_dockerfile --tag julianqz/wu_presto:main_0.1.1 \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_presto

docker build --file wu_cimm/wu_cimm_dockerfile --tag julianqz/wu_cimm:main_0.1.1 \
	--build-arg BASE_CONTAINER="wu_base:main_0.1.0" ./wu_cimm

# imgt refs

./wu_ref/wu_ref_get_imgt.sh \
	-a "/Users/jqz/Dropbox/common/germline_refs/imgt_zip/" \
	-b "/Users/jqz/Dropbox/common/germline_refs/imgt_select/" \
	-c "/Users/jqz/Dropbox/common/germline_refs/imgt_download.log" \
	-d true

Rscript ./wu_ref/wu_ref_imgt_dedup.R "202113-2" #*

# note context

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_presto:ref_0.1.1 --build-arg BASE_CONTAINER="wu_presto:main_0.1.1" \
	"/Users/jqz/Dropbox/"

docker build --file "/Users/jqz/Dropbox/wustl/code/docker/wu_ref/wu_ref_dockerfile" \
	--tag julianqz/wu_cimm:ref_0.1.1 --build-arg BASE_CONTAINER="wu_cimm:main_0.1.1" \
	"/Users/jqz/Dropbox/"

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
	--tag julianqz/ubuntu:lsf --build-arg BASE_CONTAINER="ubuntu:main" .

docker push julianqz/ubuntu:lsf

