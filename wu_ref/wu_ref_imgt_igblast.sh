#!/usr/bin/env bash
# Convert IMGT germlines sequences to IgBLAST database
#
# Author:  Julian Q Zhou
# Date:    2021-04-18
#
# Adapted from "imgt2igblast.sh" by JAVH from immcantation/scripts


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -i  Input directory containing germlines in the form:"
    echo -e "      <species>/<chain>/*.fasta."
    echo -e "  -v  Version to be noted in the name of built database."
    echo -e "      E.g.: '202113-2'."
    echo -e "  -o  Output directory for the built database."
    echo -e "  -b  Path to bin/ of IgBLAST."
    echo -e "  -h  This message."
}

# Get commandline arguments
while getopts "i:v:o:b:h" OPT; do
    case "${OPT}" in
    i)  PATH_GERM=$(realpath "${OPTARG}")
        ;;
    v)  VER="${OPTARG}"
        ;;
    o)  PATH_DB=$(realpath "${OPTARG}")
        ;;
    b)  PATH_BIN=$(realpath "${OPTARG}")
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

# so that `makeblastdb` outputs files in the database directory

cd "${PATH_DB}"

# within each species and within each chain, by segment

for SPECIES in Homo_sapiens Mus_musculus; do 

    # TODO: generalize to any species present in ${GERMDIR}

    if [[ ${SPECIES} == "Homo_sapiens" ]]; then
        SP="hs"
    else
        SP="mm"
    fi

    for CHAIN in IG TR; do

        PATH_IMGT="${PATH_GERM}/${SPECIES}/${CHAIN}/"

        for SEGMENT in V D J; do

            # concat all fasta's belonging to the same segment type (all V; all D; all J)

            NAME_CAT="concat_no_dup_${SEGMENT}.fasta"

            cat "${PATH_IMGT}"${CHAIN}?${SEGMENT}_no_dup.fasta > "${PATH_IMGT}${NAME_CAT}"
            
            # process
            # note version of IMGT reference in output files in igblast database/

            NAME_DB="imgt_no_dup_${SP}_${CHAIN}_${SEGMENT}_${VER}"

            "${PATH_BIN}/edit_imgt_file.pl" "${PATH_IMGT}${NAME_CAT}" > "${PATH_DB}/${NAME_DB}"

            "${PATH_BIN}/makeblastdb" -parse_seqids -dbtype nucl -in "${PATH_DB}/${NAME_DB}"
            
        done
    done
done

echo "Finished."
