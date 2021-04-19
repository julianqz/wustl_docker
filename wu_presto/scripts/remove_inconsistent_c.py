#!/usr/bin/env python3
"""
Remove sequences with inconsistent PRCONS and CREGION
"""
import re, sys, os
from os import path
from sys import argv

# command line input: input.fastq

# header: 
# @TAGTATTAATCTACGCA|CONSCOUNT=1,1|PRCONS=Human-IGHG|SEQORIENT=F|CREGION=Human-IGHG-InternalC

HEADER_MID_FW = "PRCONS"
HEADER_MID_RV = "CREGION"

PATTERN_MID_FW = re.compile(HEADER_MID_FW+"=[\w_-]+\|?")
PATTERN_MID_RV = re.compile(HEADER_MID_RV+"=[\w_-]+\|?")

# function to separate, if any, inconsistent pairs from consistent ones
# outputs 3 lists:
# - consistentLines: fastq lines corresponding to sequences with consistent pairs
# - inconsistentLines: fastq lines corresponding to sequences with inconsistent pairs
# - inconsistentRecords: "MID_FW, MID_RV" corresponding to sequences with inconsistent pairs
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
            
            # parse for MID_FW
            MIDfwMatch = PATTERN_MID_FW.findall(curHeader)

            if len(MIDfwMatch)==1:
                # remove MID_FW=
                MIDfw = re.sub(HEADER_MID_FW+"=", "", MIDfwMatch[0])
                # if there's |, remove it
                MIDfw = re.sub("\|", "", MIDfw)
            elif len(MIDfwMatch)==0:
                sys.exit("No match to PATTERN_MID_FW found in header")
            else:
                sys.exit("More than exactly one matches to PATTERN_MID_FW found in header")
            
            # parse for MID_RV
            MIDrvMatch = PATTERN_MID_RV.findall(curHeader)

            if len(MIDrvMatch)==1:
                # remove MID_RV=
                MIDrv = re.sub(HEADER_MID_RV+"=", "", MIDrvMatch[0])
                # if there's |, remove it
                MIDrv = re.sub("\|", "", MIDrv)
            elif len(MIDrvMatch)==0:
                sys.exit("No match to PATTERN_MID_RV found in header")
            else:
                sys.exit("More than exactly one matches to PATTERN_MID_RV found in header")

            # if not the last header
            # copy over lines containing current sequence
            if i!=(len(idxHeaders)-1):
                nextIdx = idxHeaders[i+1]
            else:
            # if the last header
            # copy over remaining lines containing sequence
                nextIdx = len(lines)

            # if matched 
            if MIDfw in MIDrv:
                consistentLines.append(curHeader)

                for seqIdx in range(idx+1, nextIdx):
                    consistentLines.append(lines[seqIdx])
            else:
                inconsistentLines.append(curHeader)
                inconsistentRecords.append(MIDfw+", "+MIDrv)

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

        print("> Output for sequences with consistent MID_FW and MID_RV: "+outputNameConsistent+"\n")
    else:
        print("> No sequence with consistent MID_FW and MID_RV found.\n")

    if numInconsistentLines>=1:
        outputNameInconsistent = outputNameStem+"_inconsistent.fastq"
        
        with open(outputNameInconsistent, 'w') as out_handle:
            out_handle.writelines(line+"\n" for line in inconsistentLines)

        print("> Output for sequences with inconsistent MID_FW and MID_RV: "+outputNameInconsistent+"\n")
    else:
        print("> No sequence with inconsistent MID_FW and MID_RV found.\n")
        
    if len(inconsistentRecords)>=1:
        outputNameRecord = outputNameStem+"_inconsistent_count.txt"
        
        # tabulate counts for MID_FW-MID_RV pairs
        recordCount = [item+", {0}".format(inconsistentRecords.count(item)) for item in set(inconsistentRecords)]
        
        with open(outputNameRecord, 'w') as out_handle:
            out_handle.writelines("MID_FW, MID_RV, COUNT\n")
            out_handle.writelines(line+"\n" for line in recordCount)

        print("> Output for counts of pairs of inconsistent MID_FW and MID_RV: "+outputNameRecord+"\n")

    print("# consistent seqs: {0} ; # inconsistent seqs: {1}.\n".format(numConsistentLines, numInconsistentLines))


# run
for f in sys.argv[1:]:

	print("\n* Starting to process {0} ...\n".format(f))

	consistentLines, inconsistentLines, inconsistentRecords = spotInconsistency(f)

	outputLines(consistentLines, inconsistentLines, inconsistentRecords, path.splitext(f)[0])

	print("* Finished processing {0}.\n".format(f))
