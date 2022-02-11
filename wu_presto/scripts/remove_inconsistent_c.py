#!/usr/bin/env python3
"""
Remove sequences with inconsistent PRCONS and CREGION
"""
import re, sys, os
from os import path
from sys import argv

# Julian Q. Zhou
# https://github.com/julianqz

# command line input: input.fastq

# header: 
# @TAGTATTAATCTACGCA|CONSCOUNT=1,1|PRCONS=Human-IGHG|SEQORIENT=F|CREGION=Human-IGHG-InternalC

HEADER_1 = "PRCONS"  #*
HEADER_2 = "CREGION" #*

PATTERN_1 = re.compile(HEADER_1+"=[\w_-]+\|?")
PATTERN_2 = re.compile(HEADER_2+"=[\w_-]+\|?")

# function to separate, if any, inconsistent pairs from consistent ones
# outputs 3 lists:
# - consistentLines: fastq lines corresponding to sequences with consistent pairs
# - inconsistentLines: fastq lines corresponding to sequences with inconsistent pairs
# - inconsistentRecords: "${HEADER_1}, ${HEADER_2}" corresponding to sequences with inconsistent pairs
def spotInconsistency(inputFasta):
    
    with open(inputFasta, "r") as in_handle:

        lines = [line.rstrip() for line in in_handle]
        consistentLines = []
        inconsistentLines = []
        inconsistentRecords = []

        # find index of lines containing "CONSCOUNT"
        # for FASTQ, can't look for "@" as header start because quality score could contain @
        idxHeaders = [idx for idx in range(len(lines)) if "CONSCOUNT" in lines[idx]]

        for i in range(len(idxHeaders)):
            idx = idxHeaders[i]

            # header
            curHeader = lines[idx]
            
            # parse for ${HEADER_1}
            match_1 = PATTERN_1.findall(curHeader)

            if len(match_1)==1:
                # remove ${HEADER_1}=
                value_1 = re.sub(HEADER_1+"=", "", match_1[0])
                # if there's |, remove it
                value_1 = re.sub("\|", "", value_1)
            elif len(match_1)==0:
                sys.exit("No match to "+HEADER_1+" found in header")
            else:
                sys.exit("More than exactly one matches to "+HEADER_1+" found in header")
            
            # parse for ${HEADER_2}
            match_2 = PATTERN_2.findall(curHeader)

            if len(match_2)==1:
                # remove ${HEADER_2}=
                value_2 = re.sub(HEADER_2+"=", "", match_2[0])
                # if there's |, remove it
                value_2 = re.sub("\|", "", value_2)
            elif len(match_2)==0:
                sys.exit("No match to "+HEADER_2+" found in header")
            else:
                sys.exit("More than exactly one matches to "+HEADER_2+" found in header")

            # if not the last header
            # copy over lines containing current sequence
            if i!=(len(idxHeaders)-1):
                nextIdx = idxHeaders[i+1]
            else:
            # if the last header
            # copy over remaining lines containing sequence
                nextIdx = len(lines)

            # if matched
            #* assumes PRCONS being a substring of CREGION if consistent
            if value_1 in value_2:
                consistentLines.append(curHeader)

                for seqIdx in range(idx+1, nextIdx):
                    consistentLines.append(lines[seqIdx])
            else:
                inconsistentLines.append(curHeader)
                inconsistentRecords.append(value_1+", "+value_2)

                for seqIdx in range(idx+1, nextIdx):
                    inconsistentLines.append(lines[seqIdx])                
                    
    return consistentLines, inconsistentLines, inconsistentRecords

# function to output files
# inputs: 3 lists from spotInconsistency(), and outputNameStem (to be used as stem for output filenames)
# outputs: files with names 
# - outputNameStem_consistent.fastq
# - outputNameStem_inconsistent.fastq
# - outputNameStem_inconsistent_count.txt
# a file will only be outputted if sequences that meet criteria of that file were found
def outputLines(consistentLines, inconsistentLines, inconsistentRecords, outputNameStem):

    numConsistentLines = sum(["CONSCOUNT" in item for item in consistentLines])
    numInconsistentLines = sum(["CONSCOUNT" in item for item in inconsistentLines])
    
    if numConsistentLines>=1:
        outputNameConsistent = outputNameStem+"_consistent.fastq"
        
        with open(outputNameConsistent, 'w') as out_handle:
            out_handle.writelines(line+"\n" for line in consistentLines)

        print("> Output for sequences with consistent "+HEADER_1+" and "+HEADER_2+": "+outputNameConsistent+"\n")
    else:
        print("> No sequence with consistent "+HEADER_1+" and "+HEADER_2+" found.\n")

    if numInconsistentLines>=1:
        outputNameInconsistent = outputNameStem+"_inconsistent.fastq"
        
        with open(outputNameInconsistent, 'w') as out_handle:
            out_handle.writelines(line+"\n" for line in inconsistentLines)

        print("> Output for sequences with inconsistent "+HEADER_1+" and "+HEADER_2+": "+outputNameInconsistent+"\n")
    else:
        print("> No sequence with inconsistent "+HEADER_1+" and "+HEADER_2+" found.\n")
        
    if len(inconsistentRecords)>=1:
        outputNameRecord = outputNameStem+"_inconsistent_count.txt"
        
        # tabulate counts for ${HEADER_1}-${HEADER_2} pairs
        recordCount = [item+", {0}".format(inconsistentRecords.count(item)) for item in set(inconsistentRecords)]
        
        with open(outputNameRecord, 'w') as out_handle:
            out_handle.writelines(HEADER_1+", "+HEADER_2+", COUNT\n")
            out_handle.writelines(line+"\n" for line in recordCount)

        print("> Output for counts of pairs of inconsistent "+HEADER_1+" and "+HEADER_2+": "+outputNameRecord+"\n")

    print("# consistent seqs: {0} ; # inconsistent seqs: {1}.\n".format(numConsistentLines, numInconsistentLines))


# run
for f in sys.argv[1:]:

	print("\n* Starting to process {0} ...\n".format(f))

	consistentLines, inconsistentLines, inconsistentRecords = spotInconsistency(f)

	outputLines(consistentLines, inconsistentLines, inconsistentRecords, path.splitext(f)[0])

	print("* Finished processing {0}.\n".format(f))
