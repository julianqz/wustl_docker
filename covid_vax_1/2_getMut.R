# same versions as used for Nature 2020
# can't install shazam 0.1.11 in R 4.0.2 b/c SDMTools not avail
#devtools::install_version("alakazam", version="0.3.0")
#devtools::install_version("shazam", version="0.1.11")

# used alakazam 1.0.2 and shazam 1.0.2 instead

library(shazam)

library(doMC)
library(foreach)
NPROC = 7
registerDoMC(NPROC) 

pathRoot = "~/Dropbox (recherche)/wustl/projects/covid_vax_1/"
pathWork = paste0(pathRoot, "JZ/mutation/")

fn = "seqDb.RData"
setwd(pathWork)
load(fn) # db

### calculate mutation frequency

EXCL_FIRST_NUM = 18

cat("\nTrimming first", EXCL_FIRST_NUM, "positions...\n")

# will return "" if EXCL_FIRST_NUM+1 > nchar (i.e. input shorter than EXCL_FIRST_NUM chars)
trimmed_SEQUENCE_IMGT = substring(text=db[["SEQUENCE_IMGT"]], first=(EXCL_FIRST_NUM+1))
trimmed_GERMLINE_IMGT = substring(text=db[["GERMLINE_IMGT"]], first=(EXCL_FIRST_NUM+1))

# regionDefinition needs to change too
IMGT_V_BY_SEGMENTS_excl = IMGT_V_BY_SEGMENTS
IMGT_V_BY_SEGMENTS_excl@boundaries = IMGT_V_BY_SEGMENTS_excl@boundaries[-(1:EXCL_FIRST_NUM)]
IMGT_V_BY_SEGMENTS_excl@seqLength = length(IMGT_V_BY_SEGMENTS_excl@boundaries)

#### nuc ####

cat("\nPost-trim, counting nucleotide mutations...\n")

mutObj_excl_nuc = foreach(i=1:nrow(db)) %dopar% {
    # obj is a list, with components "pos" and "nonN"
    # nonN is a whole number
    # pos is either a data.frame, or NA
    #     if data.frame, pos has columns "position", "R", "S", "region" 
    obj = calcObservedMutations(inputSeq=trimmed_SEQUENCE_IMGT[i], #*
                                germlineSeq=trimmed_GERMLINE_IMGT[i], #*
                                regionDefinition=IMGT_V_BY_SEGMENTS_excl, #*
                                returnRaw=T)
    return(obj)
}

# parse thru mutObj_excl_nuc

cat("\nPost-trim, parsing result from counting nucleotide mutations...\n")

mutObj_excl_nuc_parsed_lst = foreach(i=1:length(mutObj_excl_nuc)) %dopar% {
    
    o = mutObj_excl_nuc[[i]]
    
    # do not use is.na(o$pos) -- this will return multiple booleans if o$pos is data.frame
    if (is.data.frame(o$pos)) {
        nuc_R = sum(o$pos[["r"]]) # lowercase in newer shazam version
        nuc_S = sum(o$pos[["s"]]) # lowercase in newer shazam version
        nuc_RS = nuc_R+nuc_S
    } else {
        nuc_R = 0
        nuc_S = 0
        nuc_RS = 0
    }
    
    nuc_denom = o$nonN
    # otherwise the name in the returned vector ends up being nuc_denom.V
    names(nuc_denom)=NULL
    
    if (nuc_denom!=0) {
        nuc_RS_freq = nuc_RS / nuc_denom
    } else {
        nuc_RS_freq = NA
    }
    
    return( c( nuc_denom=nuc_denom, nuc_R=nuc_R, nuc_S=nuc_S, nuc_RS=nuc_RS, nuc_RS_freq=nuc_RS_freq ) )
}

#   nuc_denom nuc_R nuc_S nuc_RS nuc_RS_freq
# 1       297     8     3     11 0.037037037
mutObj_excl_nuc_parsed = data.frame( do.call(rbind, mutObj_excl_nuc_parsed_lst) )

# add trim_ to differentiate name
colnames(mutObj_excl_nuc_parsed) = paste0("trim_", colnames(mutObj_excl_nuc_parsed))

summary(mutObj_excl_nuc_parsed[["trim_nuc_RS_freq"]])

# export
db = cbind(db, mutObj_excl_nuc_parsed)

fn = "seqDbMut.RData"
setwd(pathWork)
save(db, file=fn)

db_tsv = db[, c("SEQUENCE_ID", "SAMPLE", "mAb", "CLONAL", "trim_nuc_RS_freq")]
colnames(db_tsv)[which(colnames(db_tsv)=="trim_nuc_RS_freq")] = "MF"

head(db_tsv)

fn_tsv = "mufreq_covid.tsv"
write.table(x=db_tsv, file=fn_tsv, quote=F, sep="\t",
            row.names=F, col.names=T)
