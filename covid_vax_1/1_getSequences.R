pathRoot = "~/Dropbox (recherche)/wustl/projects/covid_vax_1/"
pathNB = paste0(pathRoot, "NB/")
pathIgblast = paste0(pathRoot, "JZ/igblast/")
pathWork = paste0(pathRoot, "JZ/mutation/")

#### get seq IDs ####

library(Seurat)

setwd(pathNB)
fn = "integratedSeuratObject_subsetremoved.rds"
obj = readRDS(fn)

colnames(obj@meta.data)
head(obj@meta.data)

#  Yes  <NA> 
# 1729 10827
table(obj@meta.data[["clone_cluster"]], useNA="ifany")

#        S1   S2   S3
# Yes   338 1384    7
# <NA> 4241 3072 3514
table(obj@meta.data[["clone_cluster"]], obj@meta.data[["orig.ident"]], useNA="ifany")

# sequences of interest
idx = which(obj@meta.data[["clone_cluster"]]=="Yes")
# 1729
length(idx)

#  S1   S2   S3 
# 338 1384    7
table(obj@meta.data[["orig.ident"]][idx], useNA="ifany")

# TRUE
all.equal(rownames(obj@meta.data)[idx], obj@meta.data[["barcode"]][idx])

# 0
sum(is.na(obj@meta.data[["sequence"]][idx]))

# export
targetInfoCols = c("SEQUENCE_ID", "SAMPLE", "mAb_ID")
targetInfo = data.frame(matrix(NA, nrow=length(idx), ncol=length(targetInfoCols)))
colnames(targetInfo) = targetInfoCols

targetInfo[["SEQUENCE_ID"]] = obj@meta.data[["barcode"]][idx]
targetInfo[["SAMPLE"]] = obj@meta.data[["orig.ident"]][idx]
targetInfo[["mAb_ID"]] = obj@meta.data[["sequence"]][idx]

# all characters
sapply(1:ncol(targetInfo), function(i){class(targetInfo[, i])})

head(targetInfo)

setwd(pathWork)
#fn = "seqID.tsv"
#write.table(x=targetInfo, file=fn, quote=F, row.names=F, col.names=T)
fn = "seqID.RData"
save(targetInfo, file=fn)

rm(obj)

#### extract seqs ####

setwd(pathIgblast)

samples = c("CovidmAbsHL_corrected", "S1_filtered_contig", 
            "S2_filtered_contig", "S3_filtered_contig")

lst_tabs = vector(mode="list", length=length(samples))
names(lst_tabs) = samples

for (s in samples) {
    fn = paste0(s, "_db-pass.tab")
    tmpDb = read.table(fn, sep="\t", header=T, stringsAsFactors=F)
    
    if (s=="CovidmAbsHL_corrected") {
        tmpDb[["mAb"]] = T
        
        # heavy only
        tmpDb = tmpDb[grepl(pattern="\\.H", x=tmpDb[["SEQUENCE_ID"]]), ]
        
    } else {
        tmpDb[["mAb"]] = F
        
        # sample ID
        curSamp = substr(s, 1, 2)
        tmpDb[["SAMPLE"]] = curSamp
        
        # modify SEQUENCE_ID
        # trim "_contig_*" from SEQUENCE_ID; then add S[123]_
        # e.g. AAGCCGCTCATCTGCC-1_contig_1
        tmpDb[["SEQUENCE_ID"]] = sapply(tmpDb[["SEQUENCE_ID"]], 
                                                function(s) {
                                                    paste0(curSamp, "_", strsplit(s, "_contig")[[1]][1])
                                                })
        # heavy only
        tmpDb = tmpDb[substr(tmpDb[["V_CALL"]], 1, 4)=="IGHV", ]
    }
    
    lst_tabs[[s]] = tmpDb
    rm(tmpDb)
    
}

# 05 88 98 
#  7 19 22
table(substr(lst_tabs[["CovidmAbsHL_corrected"]][["SEQUENCE_ID"]], 1, 2), useNA="ifany")

### mAb db
mAbSampleMatch = c("05"="S3", "88"="S2", "98"="S1")

lst_tabs[["CovidmAbsHL_corrected"]][["SAMPLE"]] = mAbSampleMatch[substr(lst_tabs[["CovidmAbsHL_corrected"]][["SEQUENCE_ID"]], 1, 2)]

# S1 S2 S3 
# 22 19  7
table(lst_tabs[["CovidmAbsHL_corrected"]][["SAMPLE"]], useNA="ifany")

### combine
dbAll = do.call(rbind, lst_tabs)

#    FALSE TRUE
# S1  4861   22
# S2  4090   19
# S3  3539    7
table(dbAll[["SAMPLE"]], dbAll[["mAb"]], useNA="ifany")

### subset to seqs of interest

# wrt dbAll
idxDbAll = match(targetInfo[["SEQUENCE_ID"]], dbAll[["SEQUENCE_ID"]]) 
stopifnot(all.equal( dbAll[["SEQUENCE_ID"]][idxDbAll], 
                     targetInfo[["SEQUENCE_ID"]] ))

# each match should be the only match
for (i in 1:nrow(targetInfo)) {
    curIdx = match(targetInfo[["SEQUENCE_ID"]][i], dbAll[["SEQUENCE_ID"]])
    if (length(curIdx)!=1) {
        cat(i, targetInfo[["SEQUENCE_ID"]][i], "\n")
    }
}

db = dbAll[idxDbAll, ]

# clonally related mAbs
db[["CLONAL"]] = targetInfo[["mAb_ID"]]

# add mAb HCs
db = rbind(db, 
           cbind(dbAll[dbAll[["mAb"]], ], 
                       CLONAL=sapply(dbAll[["SEQUENCE_ID"]][dbAll[["mAb"]]],
                                     function(s){ strsplit(s, "\\.H")[[1]][1] }) ) )

# sanity check
# IGHV 
# 1777
table(substr(db[["V_CALL"]], 1, 4), useNA="ifany")

table(db[["CLONAL"]])

# export
setwd(pathWork)
fn = "seqDb.RData"
save(db, file=fn)

