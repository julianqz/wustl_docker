#!/usr/bin/env bash
# Super script to run the pRESTO pipeline on AbVitro/NEB AbSeq data
# With UMI correction as described in Jiang et al., PNAS, 2020

# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2018.09.15
# Author:  Roy Jiang
# Date:    2018.06.18
#
# Adapted by Julian Q Zhou, 2021-04-26
#
#


# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -a  Path to sample_list_${PROJ_ID}.txt"                                                        #*JZ
    echo -e "  -b  Path to directory containing all input FASTQs."                                            #*JZ
    echo -e "  -c  Suffix of Read 1 FASTQ. Sequence beginning with the C-region or J-segment).\n" \
            "      E.g. for [sample_id]_R1.fastq, the suffix should be '_R1.fastq'."                          #*JZ
    echo -e "  -d  Suffix of Read 2 FASTQ. Sequence beginning with the leader or V-segment).\n" \
            "      E.g. for [sample_id]_R2.fastq, the suffix should be '_R2.fastq'."                          #*JZ
    echo -e "  -e  Total number of sequences to subsample across samples for choosing clustering thresholds." #*JZ
    echo -e "  -f  Read 1 FASTA primer sequences.\n" \
            "      Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta."
    echo -e "  -g  Read 2 FASTA primer or template switch sequences.\n" \
            "      Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta."
    echo -e "  -i  C-region FASTA sequences for the C-region internal to the primer.\n" \
            "      If unspecified internal C-region alignment is not performed."
    echo -e "  -j  V-segment reference file.\n" \
            "      Defaults to /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta."
    echo -e "  -k  YAML file providing description fields for report generation."
    echo -e "  -l  Output directory. Will be created if it does not exist.\n" \
            "      Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -m  The mate-pair coordinate format of the raw data.\n" \
            "      Defaults to illumina."
    echo -e "  -n  Number of subprocesses for multiprocessing tools.\n" \
            "      Defaults to the available cores."
    echo -e "  -o  Path to Python script removing inconsistent C primer and internal C alignments." #*JZ
    echo -e "  -p  Path to Python script converting FASTQ to FASTA."            #*JZ
    echo -e "  -q  Boolean. Passed to ${CS_KEEP}."                              #*JZ
    echo -e "  -r  Boolean. Whether to run the pre-indexing part of pipeline."  #*JZ
    echo -e "  -s  Boolean. Whether to run the indexing part of pipeline."      #*JZ
    echo -e "  -t  Boolean. Whether to run the post-indexing part of pipeline." #*JZ
    echo -e "  -h  This message."
}

# Argument validation variables
PATH_LIST_SET=false #*JZ
N_SUBSAMPLE_SET=false  #*JZ
R1_PRIMERS_SET=false
R2_PRIMERS_SET=false
CREGION_SEQ_SET=false
VREF_SEQ_SET=false
YAML_SET=FALSE
OUTDIR_OVERALL_SET=false #*JZ
NPROC_SET=false
COORD_SET=false

# Get commandline arguments
while getopts "a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:h" OPT; do #*JZ
    case "$OPT" in
    a)  PATH_LIST=$(realpath "${OPTARG}")  #*JZ
        PATH_LIST_SET=true                 #*JZ
        ;;
    b)  PATH_INPUT=$(realpath "${OPTARG}") #*JZ
        ;;
    c)  SUFFIX_1=$OPTARG                   #*JZ
        ;; 
    d)  SUFFIX_2=$OPTARG                   #*JZ
        ;;
    e)  N_SUBSAMPLE=$OPTARG                #*JZ
        N_SUBSAMPLE_SET=true               #*JZ
        ;;
    f)  R1_PRIMERS=$OPTARG
        R1_PRIMERS_SET=true
        ;;
    g)  R2_PRIMERS=$OPTARG
        R2_PRIMERS_SET=true
        ;;
    i)  CREGION_SEQ=$OPTARG
        CREGION_SEQ_SET=true
        ;;
    j)  VREF_SEQ=$OPTARG
        VREF_SEQ_SET=true
        ;;
    k)  YAML=$OPTARG
        YAML_SET=true
        ;;
    l)  OUTDIR_OVERALL=$OPTARG  #*JZ
        OUTDIR_OVERALL_SET=true #*JZ
        ;;
    m)  COORD=$OPTARG
        COORD_SET=true
        ;;
    n)  NPROC=$OPTARG
        NPROC_SET=true
        ;;
    o)  PATH_SCRIPT_INCONSISTENT=$(realpath "${OPTARG}") #*JZ
        ;;
    p)  PATH_SCRIPT_Q2A=$(realpath "${OPTARG}")          #*JZ
        ;;
    q)  CS_KEEP=${OPTARG}                                #*JZ
        ;;                                               
    r)  BOOL_PRE=${OPTARG}                               #*JZ
        ;;
    s)  BOOL_MID=${OPTARG}                               #*JZ
        ;;
    t)  BOOL_POST=${OPTARG}                              #*JZ
        ;;
    h)  print_usage
        exit
        ;;
    \?) echo -e "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo -e "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if required arguments are not provided
if ! ${PATH_LIST_SET} || ! ${N_SUBSAMPLE_SET}; then
    echo -e "You must specify PATH_LIST_SET and N_SUBSAMPLE using the -a and -e options." >&2
    exit 1
fi

if ! ${YAML_SET}; then
    echo -e "You must specify the description file in YAML format using the -y option." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${OUTDIR_OVERALL_SET}; then
    OUTDIR_OVERALL=${OUTNAME}
fi

if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

if ! ${COORD_SET}; then
    COORD="illumina"
fi

# Check output directory permissions
if [ -e ${OUTDIR_OVERALL} ]; then
    if ! [ -w ${OUTDIR_OVERALL} ]; then
        echo -e "Output directory '${OUTDIR_OVERALL}' is not writable." >&2
        exit 1
    fi
else
    PARENTDIR=$(dirname $(realpath ${OUTDIR_OVERALL}))
    if ! [ -w ${PARENTDIR} ]; then
        echo -e "Parent directory '${PARENTDIR}' of new output directory '${OUTDIR_OVERALL}' is not writable." >&2
        exit 1
    fi
fi

# Check R1 primers
if ! ${R1_PRIMERS_SET}; then
    R1_PRIMERS="/usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta"
elif [ -e ${R1_PRIMERS} ]; then
    R1_PRIMERS=$(realpath ${R1_PRIMERS})
else
    echo -e "File '${R1_PRIMERS}' not found." >&2
    exit 1
fi

# Check R2 primers
if ! ${R2_PRIMERS_SET}; then
    R2_PRIMERS="/usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta"
elif [ -e ${R2_PRIMERS} ]; then
    R2_PRIMERS=$(realpath ${R2_PRIMERS})
else
    echo -e "File '${R2_PRIMERS}' not found." >&2
    exit 1
fi

# Check reference sequences
if ! ${VREF_SEQ_SET}; then
    VREF_SEQ="/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta"
elif [ -e ${VREF_SEQ} ]; then
    VREF_SEQ=$(realpath ${VREF_SEQ})
else
    echo -e "File '${VREF_SEQ}' not found." >&2
    exit 1
fi

# Check for C-region file
if ! ${CREGION_SEQ_SET}; then
    ALIGN_CREGION=false
elif [ -e ${CREGION_SEQ} ]; then
    ALIGN_CREGION=true
    CREGION_SEQ=$(realpath ${CREGION_SEQ})
else
    echo -e "File '${CREGION_SEQ}' not found." >&2
    exit 1
fi

# Check report yaml file
if [ -e ${YAML} ]; then
    YAML=$(realpath ${YAML})
else
    echo -e "File '${YAML}' not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=false           #*JZ
DELETE_FILES=false        #*JZ
FILTER_LOWQUAL=true
ALIGN_SETS=false
MASK_LOWQUAL=false
REPORT=true               
REMOVE_INCONSISTENT=true  #*JZ 

# FilterSeq run parameters
FS_QUAL=20
FS_MASK=30

# MaskPrimers run parameters
MP_UIDLEN=17
MP_R1_MAXERR=0.2
MP_R2_MAXERR=0.5
CREGION_MAXLEN=100
CREGION_MAXERR=0.3

# AlignSets run parameters
MUSCLE_EXEC=muscle

# BuildConsensus run parameters
BC_PRCONS_FLAG=true
BC_ERR_FLAG=true
BC_QUAL=0
BC_MINCOUNT=1
BC_MAXERR=0.1
BC_PRCONS=0.6
BC_MAXGAP=0.5

# AssemblePairs-sequential run parameters
AP_MAXERR=0.3
AP_MINLEN=8
AP_ALPHA=1e-5
AP_MINIDENT=0.5
AP_EVALUE=1e-5
AP_MAXHITS=100

# CollapseSeq run parameters
#CS_KEEP=true #*JZ now set via command line argument
CS_MISS=0


# overall log
PATH_LOG="${OUTDIR_OVERALL}/log_bcr_presto-abseq_UMI_correction_$(date '+%m%d%Y_%H%M%S').log"

N_LINES=$(wc -l < "${PATH_LIST}")
echo "N_LINES: ${N_LINES}" &> "${PATH_LOG}"


################
# Pre-indexing #
################

if $BOOL_PRE; then

    echo "Pre-indexing; by sample" &>> "${PATH_LOG}"

    for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

        # read current line ==> sample ID
        OUTNAME=$(sed "${IDX}q;d" "${PATH_LIST}") 

        echo "IDX: ${IDX}; SAMPLE: ${OUTNAME}" &>> "${PATH_LOG}"


        # sample-specific paths to R1 & R2 FASTQs
        R1_READS="${PATH_INPUT}/${OUTNAME}${SUFFIX_1}"
        R2_READS="${PATH_INPUT}/${OUTNAME}${SUFFIX_2}"

        # Check R1 reads
        if [ -e ${R1_READS} ]; then
            R1_READS=$(realpath ${R1_READS})
        else
            echo -e "File '${R1_READS}' not found." >&2
            exit 1
        fi

        # Check R2 reads
        if [ -e ${R2_READS} ]; then
            R2_READS=$(realpath ${R2_READS})
        else
            echo -e "File '${R2_READS}' not found." >&2
            exit 1
        fi

        # sample-specific folder
        OUTDIR="${OUTDIR_OVERALL}/${OUTNAME}"

        mkdir -p "${OUTDIR}"; cd "${OUTDIR}"


        # Define log files for current sample
        LOGDIR="logs"
        REPORTDIR="report"
        PIPELINE_LOG="${LOGDIR}/pipeline-presto.log"
        ERROR_LOG="${LOGDIR}/pipeline-presto.err"
        mkdir -p ${LOGDIR}
        mkdir -p ${REPORTDIR}
        echo '' > $PIPELINE_LOG
        echo '' > $ERROR_LOG

        # Check for errors
        check_error() {
            if [ -s $ERROR_LOG ]; then
                echo -e "ERROR:"
                cat $ERROR_LOG | sed 's/^/    /'
                exit 1
            fi
        }

        # Start
        PRESTO_VERSION=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
        echo -e "IDENTIFIER: ${OUTNAME}"
        echo -e "DIRECTORY: ${OUTDIR}"
        echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
        echo -e "\nSTART"
        STEP=0

        # Remove low quality reads
        if $FILTER_LOWQUAL; then
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
            FilterSeq.py quality -s $R1_READS -q $FS_QUAL --nproc $NPROC \
                --outname "${OUTNAME}-R1" --outdir . --log "${LOGDIR}/quality-1.log" \
                >> $PIPELINE_LOG  2> $ERROR_LOG
            FilterSeq.py quality -s $R2_READS -q $FS_QUAL --nproc $NPROC \
                --outname "${OUTNAME}-R2" --outdir . --log "${LOGDIR}/quality-2.log"  \
                >> $PIPELINE_LOG  2> $ERROR_LOG
            MPR1_FILE="${OUTNAME}-R1_quality-pass.fastq"
            MPR2_FILE="${OUTNAME}-R2_quality-pass.fastq"
            check_error
        else
            MPR1_FILE=$R1_READS
            MPR2_FILE=$R2_READS
        fi


        # Identify primers and UID 
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
        MaskPrimers.py score -s $MPR1_FILE -p $R1_PRIMERS --mode cut \
            --start 0 --maxerror $MP_R1_MAXERR --nproc $NPROC \
            --log "${LOGDIR}/primers-1.log" --outname "${OUTNAME}-R1" --outdir . \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        MaskPrimers.py score -s $MPR2_FILE -p $R2_PRIMERS --mode cut \
            --start $MP_UIDLEN --barcode --maxerror $MP_R2_MAXERR --nproc $NPROC \
            --log "${LOGDIR}/primers-2.log" --outname "${OUTNAME}-R2" --outdir . \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error


        # Assign UIDs to read 1 sequences
        # header still has full illumina id
        # eg: @M00990:604:000000000-JFK82:1:1101:18445:1815 1:N:0:ATCACG|PRIMER=Human-IGHG|BARCODE=CTTGTTGTTTCATATAA
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
        PairSeq.py -1 "${OUTNAME}-R1_primers-pass.fastq" -2 "${OUTNAME}-R2_primers-pass.fastq" \
            --2f BARCODE --coord $COORD >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error


        # convert header to presto format
        # [outname]_convert-pass.fastq
        # this step essential, otherwise loss of part of illumina header after AssemblePairs join
        # will cause PairSeq in post-indexing to fail
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ConvertHeaders ${COORD}"
        ConvertHeaders.py illumina \
            -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" \
            --outname "${OUTNAME}-R1_primers-pass_pair-pass" --outdir . --failed \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        ConvertHeaders.py illumina \
            -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" \
            --outname "${OUTNAME}-R2_primers-pass_pair-pass" --outdir . --failed \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error

        # Concatenate the sequence pairs end-to-end for EstimateError
        # without conversion, illumina header loses the part that comes after space (" 1:N:0:ATCACG" in the eg below)
        # eg: @M00990:604:000000000-JFK82:1:1101:18445:1815|PRIMER=Human-IGHG|BARCODE=CTTGTTGTTTCATATAA
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs join"
        AssemblePairs.py join \
            -1 "${OUTNAME}-R1_primers-pass_pair-pass_convert-pass.fastq" \
            -2 "${OUTNAME}-R2_primers-pass_pair-pass_convert-pass.fastq" \
            --1f PRIMER --2f BARCODE --coord presto --nproc $NPROC \
            --outname "${OUTNAME}-INDEX" --outdir . --log "${LOGDIR}/AssemblePairsJoin.log" \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error

        # Add the sample identity to the outputs
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders add"
        ParseHeaders.py add \
            -s "${OUTNAME}-INDEX_assemble-pass.fastq" \
            -f SAMPLE -u "${OUTNAME}" \
            --outname "${OUTNAME}-INDEX" --outdir . \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error

    done

    echo "Finished pre-indexing" &>> "${PATH_LOG}"

else
    if $FILTER_LOWQUAL; then
        STEP=5
    else
        STEP=4
    fi
fi



############
# Indexing #
############

if $BOOL_MID; then

    echo "Indexing correction; across samples" &>> "${PATH_LOG}"

    OUTNAME="JOIN"
    OUTDIR="${OUTDIR_OVERALL}/indexing/"

    mkdir -p "${OUTDIR}"; cd "${OUTDIR}"

    # Define log files
    LOGDIR="logs"
    REPORTDIR="report"
    PIPELINE_LOG="${LOGDIR}/pipeline-indexing_3.log" #&
    ERROR_LOG="${LOGDIR}/pipeline-indexing_3.err" #&
    mkdir -p ${LOGDIR}
    mkdir -p ${REPORTDIR}
    date > $PIPELINE_LOG
    date > $ERROR_LOG

    # Check for errors
    check_error() {
        if [ -s $ERROR_LOG ]; then
            echo -e "ERROR:"
            cat $ERROR_LOG | sed 's/^/    /'
            exit 1
        fi
    }

    # Start
    PRESTO_VERSION=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
    echo -e "IDENTIFIER: ${OUTNAME}"
    echo -e "DIRECTORY: ${OUTDIR}"
    echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
    echo -e "\nSTART"
    STEP_IDX=4 #&

    # concat
    # printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "cat"
    
    # echo "concat:"
    # ls "${OUTDIR_OVERALL}"/*/*-INDEX_reheader.fastq
    # cat "${OUTDIR_OVERALL}"/*/*-INDEX_reheader.fastq > JOIN.fastq

    # echo "concat:"
    # ls "${OUTDIR_OVERALL}"/*/*-R1_primers-pass_pair-pass_convert-pass.fastq
    # cat "${OUTDIR_OVERALL}"/*/*-R1_primers-pass_pair-pass_convert-pass.fastq > JOIN-R1.fastq
    
    # echo "concat:"
    # ls "${OUTDIR_OVERALL}"/*/*-R2_primers-pass_pair-pass_convert-pass.fastq
    # cat "${OUTDIR_OVERALL}"/*/*-R2_primers-pass_pair-pass_convert-pass.fastq > JOIN-R2.fastq
    

    # # subsample
    # printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "SplitSeq sample"
    # SplitSeq.py sample \
    #     -s "JOIN.fastq" -n "${N_SUBSAMPLE}" \
    #     --outname "${OUTNAME}" --outdir . \
    #     >> $PIPELINE_LOG 2> $ERROR_LOG
    # check_error

    # # 
    # printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "EstimateError barcode"
    # EstimateError.py barcode \
    #     -s "${OUTNAME}_sample1-n${N_SUBSAMPLE}.fastq" \
    #     -f BARCODE \
    #     --outname "${OUTNAME}" --outdir . \
    #     >> $PIPELINE_LOG 2> $ERROR_LOG
    # check_error


    # TABLE_UID="${OUTDIR}/${OUTNAME}_threshold-barcode.tab"
    
    # CMD_UID="import pandas as pd; threshold=1-pd.read_table('${TABLE_UID}', index_col='TYPE')['THRESHOLD']['ALL']; print(0.8 if threshold<0.8 else threshold)"

    # UID_THRESHOLD_PERCENT=$(echo -e "${CMD_UID}" | python3)

    # echo "thresh-barcode: ${UID_THRESHOLD_PERCENT}" &>> "${PATH_LOG}"


    # # cluster barcodes
    # printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "ClusterSets barcode"
    # ClusterSets.py barcode \
    #     -s "JOIN.fastq" \
    #     -f BARCODE \
    #     -k INDEX_UID \
    #     --ident "${UID_THRESHOLD_PERCENT}" \
    #     --cluster cd-hit-est \
    #     --prefix "u" \
    #     --nproc "${NPROC}" \
    #     --outname "${OUTNAME}-uid" --outdir . \
    #     >> $PIPELINE_LOG 2> $ERROR_LOG
    # check_error

    # subsample
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "SplitSeq sample"
    SplitSeq.py sample \
        -s "${OUTNAME}-uid_cluster-pass.fastq" -n "${N_SUBSAMPLE}" \
        --outname "${OUTNAME}-uid" --outdir . \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error    
    
    # defaults:
    #   -n: min num of seqs needed to consider a set: 20
    #   --mode: method for determining consensus sequence: freq (alternative is qual)
    #   --freq: MIN_FREQ: 0.6
    #   -q: MIN_QUAL: 0
    #   --maxdiv: specify to calc nucleotide diversity of each read group and exclude
    #             groups which exceed the given diversity thresold: None

    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "EstimateError set"
    EstimateError.py set \
        -s "${OUTNAME}-uid_sample1-n${N_SUBSAMPLE}.fastq" \
        -f INDEX_UID \
        --nproc "${NPROC}" \
        --outname "${OUTNAME}" --outdir . --log "${LOGDIR}/EstimateErrorSet.log" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error


    TABLE_SET="${OUTDIR}/${OUTNAME}_threshold-set.tab"

    CMD_SET="import pandas as pd; threshold=1-pd.read_table('${TABLE_SET}', index_col='TYPE')['THRESHOLD']['ALL']; print(0.8 if threshold<0.8 else threshold)"

    SEQ_THRESHOLD_PERCENT=$(echo -e "${CMD_SET}" | python3)

    echo "thresh-set: ${SEQ_THRESHOLD_PERCENT}" &>> "${PATH_LOG}"


    # cluster sequences within the UID groups across all samples
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "ClusterSets set"
    ClusterSets.py set \
        -s "${OUTNAME}-uid_cluster-pass.fastq" \
        -f INDEX_UID \
        -k INDEX_SEQ \
        --ident "${SEQ_THRESHOLD_PERCENT}" \
        --failed \
        --cluster cd-hit-est \
        --prefix "s" \
        --nproc "${NPROC}" \
        --outname "${OUTNAME}-seq" --outdir . --log "${LOGDIR}/ClusterSetsSet.log" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error

    # collapse the INDEX_UID and INDEX_SEQ to create a final "true" barcode
    
    # default:
    #   --act {min,max,sum,set,cat}: None
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "ParseHeaders merge"
    ParseHeaders.py merge \
        -s "${OUTNAME}-seq_cluster-pass.fastq" \
        -f INDEX_UID INDEX_SEQ \
        -k INDEX_NEW \
        --failed \
        --outname "${OUTNAME}-seq-merge" --outdir . \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error

    # make a copy of SAMPLE
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "ParseHeaders copy"
    ParseHeaders.py copy \
        -s "${OUTNAME}"-seq-merge_reheader.fastq \
        -f SAMPLE -k SAMPLE_ORIG --act set \
        --failed \
        --outname "${OUTNAME}-seq-copy" --outdir . \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error

    # resolve collisions
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "UnifyHeaders consensus"
    UnifyHeaders.py consensus \
        -s "${OUTNAME}"-seq-copy_reheader.fastq \
        -f INDEX_NEW \
        -k SAMPLE \
        --failed \
        --nproc "${NPROC}" \
        --outname "${OUTNAME}" --outdir . --log "${LOGDIR}/UnifyHeaders.log" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error

    # divide into samples
    # ${OUTNAME}_SAMPLE-${sample_id}.fastq
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP_IDX)) 24 "SplitSeq group"
    SplitSeq.py group \
        -s "${OUTNAME}_unify-pass.fastq" \
        -f SAMPLE \
        --outname "${OUTNAME}" --outdir . \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error


    echo "Finished indexing correction" &>> "${PATH_LOG}"

fi



#################
# Post-indexing #
#################

if $BOOL_POST; then

    echo "Post-indexing; by sample" &>> "${PATH_LOG}"

    # indexing folder
    OUTDIR_MID="${OUTDIR_OVERALL}/indexing"

    # keep a constant copy of the starting step outside for loop
    # otherwise it will keep growing for later samples
    STEP_INIT=${STEP}


    for ((IDX=1; IDX<=${N_LINES}; IDX++)); do

        # read current line ==> sample id
        OUTNAME=$(sed "${IDX}q;d" "${PATH_LIST}") 

        echo "IDX: ${IDX}; SAMPLE: ${OUTNAME}" &>> "${PATH_LOG}"


        # sample-specific folder
        OUTDIR="${OUTDIR_OVERALL}/${OUTNAME}"

        mkdir -p "${OUTDIR}"; cd "${OUTDIR}"


        # Define log files for current sample
        LOGDIR="logs"
        REPORTDIR="report"
        PIPELINE_LOG="${LOGDIR}/pipeline-presto.log"
        ERROR_LOG="${LOGDIR}/pipeline-presto.err"
        mkdir -p ${LOGDIR}
        mkdir -p ${REPORTDIR}

        if [ -s $PIPELINE_LOG ]; then
            echo '' >> $PIPELINE_LOG
        else
            echo '' > $PIPELINE_LOG
        fi

        if [ -s $ERROR_LOG ]; then
            echo '' >> $ERROR_LOG
        else
            echo '' > $ERROR_LOG
        fi

        # Check for errors
        check_error() {
            if [ -s $ERROR_LOG ]; then
                echo -e "ERROR:"
                cat $ERROR_LOG | sed 's/^/    /'
                exit 1
            fi
        }

        # Start        
        PRESTO_VERSION=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
        echo -e "IDENTIFIER: ${OUTNAME}"
        echo -e "DIRECTORY: ${OUTDIR}"
        echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
        echo -e "\nSTART"
        #STEP=$(($STEP_INIT+1))
        STEP="${STEP_INIT}"


        # 
        # ${OUTNAME}-INDEX-R1-1_pair-pass.fastq
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
        PairSeq.py \
            -1 "${OUTDIR_MID}/JOIN-R1.fastq" \
            -2 "${OUTDIR_MID}/JOIN_SAMPLE-${OUTNAME}.fastq" \
            --2f INDEX_UID INDEX_SEQ INDEX_NEW SAMPLE_ORIG \
            --coord "presto" \
            --failed \
            --outname "${OUTNAME}-INDEX-R1" --outdir "${OUTDIR_MID}" \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error

        #rm "${OUTDIR_MID}/${OUTNAME}-INDEX-R1-2_pair-pass.fastq"
        


        # ${OUTNAME}-INDEX-R2-1_pair-pass.fastq
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
        PairSeq.py \
            -1 "${OUTDIR_MID}/JOIN-R2.fastq" \
            -2 "${OUTDIR_MID}/JOIN_SAMPLE-${OUTNAME}.fastq" \
            --2f INDEX_UID INDEX_SEQ INDEX_NEW SAMPLE_ORIG \
            --coord "presto" \
            --failed \
            --outname "${OUTNAME}-INDEX-R2" --outdir "${OUTDIR_MID}" \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error

        #rm "${OUTDIR_MID}/${OUTNAME}-INDEX-R2-2_pair-pass.fastq"
        

        cp "${OUTDIR_MID}/JOIN_SAMPLE-${OUTNAME}.fastq" \
           "${OUTDIR}/${OUTNAME}-INDEX.fastq"


        cp "${OUTDIR_MID}/${OUTNAME}-INDEX-R1-1_pair-pass.fastq" \
            "${OUTDIR}/${OUTNAME}-INDEX-R1.fastq"

        BCR1_FILE="${OUTNAME}-INDEX-R1.fastq"


        cp "${OUTDIR_MID}/${OUTNAME}-INDEX-R2-1_pair-pass.fastq" \
            "${OUTDIR}/${OUTNAME}-INDEX-R2.fastq"

        BCR2_FILE="${OUTNAME}-INDEX-R2.fastq"


        #* IMPORTANT: replace --bf BARCODE with INDEX_NEW
        #* add --cf BARCODE to keep track of BARCODE (otherwise won't propogate)

        # Build UID consensus sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
        if $BC_ERR_FLAG; then
            if $BC_PRCONS_FLAG; then
                BuildConsensus.py -s $BCR1_FILE --bf INDEX_NEW --pf PRIMER --prcons $BC_PRCONS \
                    --cf BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG --act set set set set \
                    -n $BC_MINCOUNT -q $BC_QUAL --maxerror $BC_MAXERR --maxgap $BC_MAXGAP \
                    --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
                    --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
            else
                BuildConsensus.py -s $BCR1_FILE --bf INDEX_NEW --pf PRIMER \
                    --cf BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG --act set set set set \
                    -n $BC_MINCOUNT -q $BC_QUAL --maxerror $BC_MAXERR --maxgap $BC_MAXGAP \
                    --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
                    --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
            fi

            BuildConsensus.py -s $BCR2_FILE --bf INDEX_NEW --pf PRIMER \
                --cf BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG --act set set set set \
                -n $BC_MINCOUNT -q $BC_QUAL --maxerror $BC_MAXERR --maxgap $BC_MAXGAP \
                --nproc $NPROC --log "${LOGDIR}/consensus-2.log" \
                --outname "${OUTNAME}-R2" >> $PIPELINE_LOG 2> $ERROR_LOG
        else
            if $BC_PRCONS_FLAG; then
                BuildConsensus.py -s $BCR1_FILE --bf INDEX_NEW --pf PRIMER --prcons $BC_PRCONS \
                    --cf BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG --act set set set set \
                    -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP \
                    --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
                    --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
            else
                BuildConsensus.py -s $BCR1_FILE --bf INDEX_NEW --pf PRIMER \
                    --cf BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG --act set set set set \
                    -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP \
                    --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
                    --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
            fi

            BuildConsensus.py -s $BCR2_FILE --bf INDEX_NEW --pf PRIMER \
                --cf BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG --act set set set set \
                -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP \
                --nproc $NPROC --log "${LOGDIR}/consensus-2.log" \
                --outname "${OUTNAME}-R2" >> $PIPELINE_LOG 2> $ERROR_LOG
        fi
        check_error

        #* BuildConsensus changes coord into presto (JZ)

        # Syncronize read files
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
        PairSeq.py -1 "${OUTNAME}-R1_consensus-pass.fastq" -2 "${OUTNAME}-R2_consensus-pass.fastq" \
            --coord presto >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error


        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
        if $BC_PRCONS_FLAG; then
            PRFIELD="PRCONS"
        else
            PRFIELD="PRIMER"
        fi

        #* added more fields to copy over
        AssemblePairs.py sequential -1 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
            -2 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" -r $VREF_SEQ \
            --coord presto --rc tail \
            --1f CONSCOUNT BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG \
            --2f $PRFIELD CONSCOUNT PRFREQ BARCODE INDEX_UID INDEX_SEQ SAMPLE_ORIG \
            --minlen $AP_MINLEN --maxerror $AP_MAXERR --alpha $AP_ALPHA --scanrev \
            --minident $AP_MINIDENT --evalue $AP_EVALUE --maxhits $AP_MAXHITS --aligner blastn \
            --failed \
            --nproc $NPROC --log "${LOGDIR}/assemble.log" \
            --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
        PH_FILE="${OUTNAME}_assemble-pass.fastq"
        check_error


        # Mask low quality positions
        if $MASK_LOWQUAL; then
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
            FilterSeq.py maskqual -s $PH_FILE -q $FS_MASK --nproc $NPROC \
                --outname "${OUTNAME}-MQ" --log "${LOGDIR}/maskqual.log" \
                >> $PIPELINE_LOG 2> $ERROR_LOG
            PH_FILE="${OUTNAME}-MQ_maskqual-pass.fastq"
            check_error
        fi


        if $ALIGN_CREGION; then
            # Annotate with internal C-region
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
            CREGION_FIELD="CREGION"
            MaskPrimers.py align -s $PH_FILE -p $CREGION_SEQ \
                --maxlen $CREGION_MAXLEN --maxerror $CREGION_MAXERR \
                --mode tag --revpr --skiprc --pf $CREGION_FIELD \
                --log "${LOGDIR}/cregion.log" --outname "${OUTNAME}-CR" --nproc $NPROC \
                >> $PIPELINE_LOG 2> $ERROR_LOG
            PH_FILE="${OUTNAME}-CR_primers-pass.fastq"
            check_error
        else
            CREGION_FIELD=""
        fi


        #*JZ next block
        if $ALIGN_CREGION; then
            if $REMOVE_INCONSISTENT; then
                
                # outputs:
                # ${OUTNAME}-CR_primers-pass_consistent.fastq
                # ${OUTNAME}-CR_primers-pass_inconsistent.fastq (if any)
                # ${OUTNAME}-CR_primers-pass_inconsistent_count.txt (if any)

                # ${PH_FILE} expected to be ${OUTNAME}-CR_primers-pass.fastq
                python3 "${PATH_SCRIPT_INCONSISTENT}" "${PH_FILE}"
                    
                COUNT_FASTQ=`grep -c 'CONSCOUNT' ${OUTNAME}-CR_primers-pass_consistent.fastq`
                echo -e "consistent PRCONS and CREGION: ${COUNT_FASTQ}" 

                if [[ (-s "${OUTNAME}-CR_primers-pass_inconsistent.fastq") ]]; then
                    COUNT_FASTQ=`grep -c 'CONSCOUNT' ${OUTNAME}-CR_primers-pass_inconsistent.fastq`
                    echo -e "inconsistent PRCONS and CREGION: ${COUNT_FASTQ}"
                fi

                PH_FILE="${OUTNAME}-CR_primers-pass_consistent.fastq"

            fi
        fi


        # Rewrite header with minimum of CONSCOUNT
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
        ParseHeaders.py collapse -s $PH_FILE -f CONSCOUNT --act min \
            --outname "${OUTNAME}-final" > /dev/null 2> $ERROR_LOG
        mv "${OUTNAME}-final_reheader.fastq" "${OUTNAME}-final_total.fastq"
        check_error


        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        if $CS_KEEP; then
            CollapseSeq.py -s "${OUTNAME}-final_total.fastq" -n $CS_MISS \
            --uf PRCONS $CREGION_FIELD --cf CONSCOUNT --act sum --inner \
            --keepmiss --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
        else
            CollapseSeq.py -s "${OUTNAME}-final_total.fastq" -n $CS_MISS \
            --uf PRCONS $CREGION_FIELD --cf CONSCOUNT --act sum --inner \
            --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
        fi
        check_error


        # Filter to sequences with at least 2 supporting sources
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
        SplitSeq.py group -s "${OUTNAME}-final_collapse-unique.fastq" -f CONSCOUNT --num 2 \
            >> $PIPELINE_LOG 2> $ERROR_LOG
        check_error

        #* can't have "-" (dash) in shell variable name
        NROW_final_total=$((`wc -l < "${OUTNAME}-final_total.fastq"`))
        NROW_final_collapse_unique=$((`wc -l < "${OUTNAME}-final_collapse-unique.fastq"`))
        NROW_final_collapse_unique_atleast_2=$((`wc -l < "${OUTNAME}-final_collapse-unique_atleast-2.fastq"`))

        # Create table of final repertoire
        #* added if loop to catch empty .fastq
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
        if [[ $NROW_final_total > 0 ]]; then
            ParseHeaders.py table -s "${OUTNAME}-final_total.fastq" \
                -f ID PRCONS $CREGION_FIELD CONSCOUNT --outname "final-total" \
                --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
        fi
        if [[ $NROW_final_collapse_unique > 0 ]]; then
            ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique.fastq" \
                -f ID PRCONS $CREGION_FIELD CONSCOUNT DUPCOUNT --outname "final-unique" \
                --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
        fi
        if [[ $NROW_final_collapse_unique_atleast_2 > 0 ]]; then
            ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique_atleast-2.fastq" \
                -f ID PRCONS $CREGION_FIELD CONSCOUNT DUPCOUNT --outname "final-unique-atleast2" \
                --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
        fi
        check_error


        # Process log files
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
        if $FILTER_LOWQUAL; then
            ParseLog.py -l "${LOGDIR}/quality-1.log" "${LOGDIR}/quality-2.log" -f ID QUALITY \
                --outdir ${LOGDIR} > /dev/null &
        fi
        ParseLog.py -l "${LOGDIR}/primers-1.log" "${LOGDIR}/primers-2.log" -f ID BARCODE PRIMER ERROR \
            --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
        ParseLog.py -l "${LOGDIR}/consensus-1.log" "${LOGDIR}/consensus-2.log" \
            -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR \
            --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
        ParseLog.py -l "${LOGDIR}/assemble.log" \
            -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
            --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
        if $MASK_LOWQUAL; then
            ParseLog.py -l "${LOGDIR}/maskqual.log" -f ID MASKED \
                --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
        fi
        if $ALIGN_CREGION; then
            ParseLog.py -l "${LOGDIR}/cregion.log" -f ID PRIMER ERROR \
                --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
        fi
        wait
        check_error

        # Generate pRESTO report
        if $REPORT; then
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Generating report"
            REPORT_SCRIPT="buildReport(\"${LOGDIR}\", sample=\"${OUTNAME}\", output_dir=\"${REPORTDIR}\", template=\"AbSeqV3\", config=\"${YAML}\", quiet=FALSE, format=\"html\")" #*JZ
            Rscript -e "library(prestor); ${REPORT_SCRIPT}" > ${REPORTDIR}/report.out 2> ${REPORTDIR}/report.err
            #* rename prestor report
            #* name generated from yaml is the same for all samples b/c one yaml file used for entire project
            #* error: `mv "${REPORTDIR}"/*.html "${REPORTDIR}/${OUTNAME}_prestor.html"` 
            #*         mv: target 'report/sample1_prestor.html' is not a directory
            cd "${REPORTDIR}"
            mv *.html "${OUTNAME}_prestor.html"
        fi

        #* convert fastq to fasta
        cd "${OUTDIR}"
        echo "Converting output fastq to fasta"
        "${PATH_SCRIPT_Q2A}" "${OUTNAME}-final_collapse-unique_atleast-2.fastq"
        "${PATH_SCRIPT_Q2A}" "${OUTNAME}-final_collapse-unique.fastq"


        # Zip or delete intermediate and log files
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
        LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
        FILTER_FILES="$(basename ${R1_READS})\|$(basename ${R2_READS})\|$(basename ${R1_PRIMERS})\|$(basename ${R2_PRIMERS})"
        FILTER_FILES+="\|final_total.fastq\|final_collapse-unique.fastq\|final_collapse-unique_atleast-2.fastq"
        TEMP_FILES=$(ls *.fastq  2>/dev/null | grep -v ${FILTER_FILES})
        if $ZIP_FILES; then
            tar -zcf log_files.tar.gz $LOG_FILES
            tar -zcf temp_files.tar.gz $TEMP_FILES
        fi
        if $DELETE_FILES; then
            rm $TEMP_FILES
            rm $LOG_FILES
        fi


        # End
        printf "DONE\n\n"

    done

    echo "Finished post-indexing" &>> "${PATH_LOG}"

fi
