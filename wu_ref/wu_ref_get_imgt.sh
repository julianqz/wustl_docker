#!/usr/bin/env bash
# Download the latest IMGT germline references and tag with release version
#
# Author: Julian Q Zhou
# Date:   2021-04-18


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -a  Directory to store downloaded zip file."
    echo -e "  -b  Directory to store unzipped files."
    echo -e "  -c  Path of download log (e.g. .../.../[name].log)."
    echo -e "  -d  Boolean. Whether to store unzipped files for human and mouse only."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "a:b:c:d:h" OPT; do
    case "${OPT}" in
    a)  PATH_SAVE_ZIP=$(realpath "${OPTARG}")
		;;
    b)  PATH_SAVE_UNZIP=$(realpath "${OPTARG}")
        ;;
    c)  PATH_LOG=$(realpath "${OPTARG}")
		;;
	d)  HM_ONLY="${OPTARG}"
		;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done


# release notes
URL_IMGT_NOTES="http://www.imgt.org/IMGT_vquest/data_releases/"

# ftp site
URL_IMGT_DOWNLOAD="http://www.imgt.org/download/V-QUEST/"

# .txt containing version and release date
NAME_TXT="IMGT_vquest_release.txt"
URL_IMGT_TXT="${URL_IMGT_DOWNLOAD}${NAME_TXT}"

# zip file
NAME_ZIP="IMGT_V-QUEST_reference_directory.zip"
URL_IMGT_ZIP="${URL_IMGT_DOWNLOAD}${NAME_ZIP}"


##### obtain latest version number from .txt

# download .txt
echo "Getting latest version number from IMGT download site."
echo "Downloading ${NAME_TXT}."

wget --quiet -O "${PATH_SAVE_ZIP}/${NAME_TXT}" "${URL_IMGT_TXT}"

# save content of .txt as a variable
# e.g. "202113-2 (29 March 2021)"
TXT=$(cat "${PATH_SAVE_ZIP}/${NAME_TXT}")

# parse .txt for just the version number
# ${string%%substring}: deletes longest match of $substring from back of $string
# e.g. "202113-2"
VER=$(echo ${TXT%% (*})
echo "Latest version: ${VER}. Do manually confirm that this is correct."

# remove .txt
echo "Removing ${NAME_TXT}."
rm "${PATH_SAVE_ZIP}/${NAME_TXT}"


# .zip to be downloaded
# e.g. IMGT_vquest_release202113-2.zip
NAME_ZIP_VER="IMGT_vquest_release${VER}.zip"

# name of folder that .zip is to be unzippped into
# e.g. IMGT_vquest_release202113-2
NAME_UNZIP_VER="IMGT_vquest_release${VER}"


# only proceed if neither NAME_ZIP_VER nor NAME_UNZIP_VER exists

if [[ ! (-s "${PATH_SAVE_ZIP}/${NAME_ZIP_VER}") && ! (-d "${PATH_SAVE_UNZIP}/${NAME_UNZIP_VER}") ]]; then
	
	##### download .zip & rename with version
	
	echo "Downloading ${URL_IMGT_ZIP} as ${NAME_ZIP_VER}."
	wget --quiet -O "${PATH_SAVE_ZIP}/${NAME_ZIP_VER}" "${URL_IMGT_ZIP}"

	##### update download log

	echo "Updating ${PATH_LOG}"

	DT=$(date '+%Y-%m-%d')
	MSG="Release: ${TXT} | Downloaded: ${DT}"

	if [[ ! (-s "${PATH_LOG}" ) ]]; then
		# if log does not exist, create
		echo "${MSG}" &> "${PATH_LOG}"
	else
		# if exists, append
		echo "${MSG}" &>> "${PATH_LOG}"
	fi

	##### unzip

	echo "Unzipping ${NAME_ZIP_VER} into ${NAME_UNZIP_VER}."
	
	# -n: never overwrite existing files
	unzip -n -q "${PATH_SAVE_ZIP}/${NAME_ZIP_VER}" -d "${PATH_SAVE_UNZIP}"
	
	# original name of unzipped folder
	NAME_UNZIP_ORIG="IMGT_V-QUEST_reference_directory"

	if $HM_ONLY; then
		# keep only human and mouse folders in IMGT_vquest_${VER}

		echo "Removing all unzipped folders except for human and mouse folders."
		
		mkdir -p "${PATH_SAVE_UNZIP}/${NAME_UNZIP_VER}"
		
		mv "${PATH_SAVE_UNZIP}/${NAME_UNZIP_ORIG}/Homo_sapiens" \
		   "${PATH_SAVE_UNZIP}/${NAME_UNZIP_ORIG}/Mus_musculus" \
		   "${PATH_SAVE_UNZIP}/${NAME_UNZIP_VER}"
		
		rm -r "${PATH_SAVE_UNZIP}/${NAME_UNZIP_ORIG}"
	else
		# keep all folders in IMGT_vquest_${VER}

		mv "${PATH_SAVE_UNZIP}/${NAME_UNZIP_ORIG}" "${PATH_SAVE_UNZIP}/${NAME_UNZIP_VER}"
	fi


	echo "Finished downloading IMGT reference release${VER}."

else
	echo "${NAME_ZIP_VER} and/or ${NAME_UNZIP_VER}/ already exists. Skipped."
fi
