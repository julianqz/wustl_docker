#!/usr/bin/env bash
# Convert IMGT germlines sequences to IgBLAST database
#
# Author:  Julian Q. Zhou
# https://github.com/julianqz
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

#for SPECIES in Homo_sapiens Mus_musculus Mus_musculus_C57BL6; do 
for SPECIES in Homo_sapiens Mus_musculus Mus_musculus_C57BL6 H2L2; do # 0.3.2b

    # TODO: generalize to any species present in ${PATH_GERM}

    if [[ ${SPECIES} == "Homo_sapiens" ]]; then
        SP="hs"
    elif [[ ${SPECIES} == "Mus_musculus" ]]; then
        SP="mm"
    elif [[ ${SPECIES} == "Mus_musculus_C57BL6" ]]; then
        SP="c57bl6"
    elif [[ ${SPECIES} == "H2L2" ]]; then
        SP="h2l2"
    fi

    PATH_SPECIES="${PATH_GERM}/${SPECIES}/"

    # if species directory present
    if [[ -d "${PATH_SPECIES}" ]]; then

        for CHAIN in IG TR; do

            PATH_IMGT="${PATH_SPECIES}${CHAIN}/"

            # if chain directory present
            if [[ -d "${PATH_IMGT}" ]] ; then

                #for SEGMENT in V D J; do
                for SEGMENT in V D J C; do # 0.3.2b

                    NAME_CAT="concat_no_dup_${SEGMENT}.fasta"

                    # note version of IMGT reference in output files in igblast database/
                    NAME_DB="imgt_no_dup_${SP}_${CHAIN}_${SEGMENT}_${VER}"

                    if [ -s "${PATH_IMGT}${NAME_CAT}" ]; then

                        # if concatenated file already exists (produced by wu_ref_imgt_dedup.R)
                        echo "Using existing ${PATH_IMGT}${NAME_CAT}"

                        "${PATH_BIN}/edit_imgt_file.pl" "${PATH_IMGT}${NAME_CAT}" > "${PATH_DB}/${NAME_DB}"

                        "${PATH_BIN}/makeblastdb" -parse_seqids -dbtype nucl -in "${PATH_DB}/${NAME_DB}"

                    else

                        # concat all fasta's belonging to the same segment type (all V; all D; all J; all C)

                        # ${CHAIN}?${SEGMENT} captures, eg. IG[HKL]V

                        # https://stackoverflow.com/questions/6363441/check-if-a-file-exists-with-a-wildcard-in-a-shell-script

                        if compgen -G "${PATH_IMGT}${CHAIN}?${SEGMENT}_no_dup.fasta" > /dev/null; then

                            echo "Performing concatenation to produce ${NAME_CAT}"

                            cat "${PATH_IMGT}"${CHAIN}?${SEGMENT}_no_dup.fasta > "${PATH_IMGT}${NAME_CAT}"
                            
                            # process
                            "${PATH_BIN}/edit_imgt_file.pl" "${PATH_IMGT}${NAME_CAT}" > "${PATH_DB}/${NAME_DB}"

                            # IMPORTANT
                            # will get fatal errors such as the ones below if omitting -parse_seqids
                            # - WORKER: T7 BATCH # 5 CEXCEPTION: Attempt to access NULL pointer.
                            # - terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
                            #     what():  basic_ios::clear: iostream error
                            
                            "${PATH_BIN}/makeblastdb" -parse_seqids -dbtype nucl -in "${PATH_DB}/${NAME_DB}"
         
                        else
                            echo "No ${SEGMENT}_no_dup.fasta exists for ${CHAIN}; skipped"
                        fi
                    fi

                done
            fi
        done
    fi
done

echo "Finished."
