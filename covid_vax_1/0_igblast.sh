#!/bin/bash

NPROC=7

# change-o v0.4.3

PATH_ROOT="~/Dropbox (recherche)/wustl/projects/covid_vax_1/" #*
PATH_DATA="${PATH_ROOT}NB/"
PATH_IGBLAST="${PATH_ROOT}JZ/igblast/"

PATH_LOG="${PATH_IGBLAST}igblast_allSamples_$(date '+%m%d%Y_%H%M%S').log"
PATH_COUNT="${PATH_IGBLAST}igblast_bySample_count_$(date '+%m%d%Y_%H%M%S').txt"

GERM_TAG="ref_201931-4_1Aug2019" #*
FN_IMGTVDJ="/Users/jkewz/Dropbox (recherche)/yale/projects/germ/ref_201931-4_1Aug2019/human_IG/" #IG*.fasta" #*

DIR_IGDATA="/Users/jkewz/apps/igblast/ncbi-igblast-1.14.0" #*
PATH_EXEC="/Users/jkewz/apps/igblast/ncbi-igblast-1.14.0/bin/igblastn" #*


#* run control
# TRUE/FALSE
RUN_IGBLAST=TRUE
RUN_MAKEDB=TRUE

# start log file
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "$dt" &> "${PATH_LOG}"
echo "germline" $GERM_TAG >> "${PATH_LOG}"
echo "bash" $BASH_VERSION >> "${PATH_LOG}"
#parallel --version | head -n 1 >> "${PATH_LOG}"
"${PATH_EXEC}" -version >> "${PATH_LOG}"
MakeDb.py --version >> "${PATH_LOG}"
python3 -V >> "${PATH_LOG}"


SAMP_ARRAY=(CovidmAbsHL_corrected S1_filtered_contig S2_filtered_contig S3_filtered_contig)

for SAMP in ${SAMP_ARRAY[@]}; do
	
	echo "####### ${SAMP}" >> "${PATH_LOG}"

	PATH_FASTA="${PATH_DATA}${SAMP}.fasta"
	
	# sample-specific log
	PATH_LOG_SAMPLE="${PATH_IGBLAST}igblast_${SAMP}_$(date '+%m%d%Y_%H%M%S').log"

	#######################
	##### run igblast #####
	#######################

	if [[ ${RUN_IGBLAST} == "TRUE" ]]; then

		if [[ (-s "${PATH_FASTA}" ) ]]; then

			# count reads in cell ranger filtered contig fasta
			COUNT_CR=`grep -c '^>' "${PATH_FASTA}"`
			echo -e "${SAMP}_CR\t${COUNT_CR}" >> "${PATH_COUNT}"

			echo "##### run IgBLAST #####" >> "${PATH_LOG_SAMPLE}"

			# output: [outname]_igblast.fmt7

			# vdb, ddb, jdb: name of the custom V/D/J reference in IgBLAST database folder
			# if unspecified, default to imgt_<organism>_<loci>_v/d/j
			AssignGenes.py igblast \
				-s "${PATH_FASTA}" \
				-b "${DIR_IGDATA}" \
				--exec "${PATH_EXEC}" \
				--organism "human" \
				--loci "ig" \
				--vdb "imgt_human_IG_V_201931-4" \
				--ddb "imgt_human_IG_D_201931-4" \
				--jdb "imgt_human_IG_J_201931-4" \
				--format "blast" \
				--outname "${SAMP}" \
				--outdir "${PATH_IGBLAST}" \
				--nproc "${NPROC}" \
				>> "${PATH_LOG_SAMPLE}"

		else
			echo "${PATH_FASTA} does not exist" >> "${PATH_LOG_SAMPLE}"
		fi

	fi

	#######################################
	##### process output from IgBLAST #####
	#######################################

	if [[ ${RUN_MAKEDB} == "TRUE" ]]; then

		PATH_IGBLAST_SAMPLE="${PATH_IGBLAST}${SAMP}_igblast.fmt7"

		if [[ (-s "${PATH_IGBLAST_SAMPLE}" ) ]]; then
			
			echo "##### process output from IgBLAST #####" >> "${PATH_LOG_SAMPLE}"

			# default: partial=False; asis-id=False
			# --partial: if specified, include incomplete V(D)J alignments in the pass file instead of the fail file

			# v0.4.6: Combined the extended field arguments of all subcommands 
			#   (--scores, --regions, --cdr3, and --junction) into a single --extended argument

			# IMPORTANT: do not put "" around FN_IMGTVDJ (otherwise * will be interpreted as is)
			# output: [outname]_db-pass.tab
			MakeDb.py igblast \
				-i "${PATH_IGBLAST_SAMPLE}" \
				-s "${PATH_FASTA}" \
				-r "${FN_IMGTVDJ}"IG*.fasta \
				--extended \
				--failed \
				--partial \
				--format changeo \
				--outname "${SAMP}" \
				--outdir "${PATH_IGBLAST}" \
				>> "${PATH_LOG_SAMPLE}"

			# count
			COUNT_LINES=$(wc -l "${PATH_IGBLAST}${SAMP}_db-pass.tab" | awk '{print $1}')
			COUNT_SEQS=$((COUNT_LINES - 1))
			echo -e "${SAMP}_MakeDb\t${COUNT_SEQS}" >> "${PATH_COUNT}"

		else
			echo "${PATH_IGBLAST_SAMPLE} does not exist" >> "${PATH_LOG_SAMPLE}"
		fi

	fi

done

echo "DONE" >> "${PATH_LOG}"
