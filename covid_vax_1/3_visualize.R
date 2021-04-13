pathRoot = "~/Dropbox (recherche)/wustl/projects/covid_vax_1/"
pathWork = paste0(pathRoot, "JZ/mutation/")

setwd(pathWork)

library(beeswarm) # v0.2.3

#### mutation data ####

fn = "mufreq_covid.tsv"
db = read.table(fn, header=T, sep="\t", stringsAsFactors=F)

fn_nat = "mufreq_nature.tsv"
db_nat = read.table(fn_nat, header=T, sep="\t", stringsAsFactors=F)

nat_naive = db_nat[["MF"]][db_nat[["SUBSET"]]=="Naive"]
# "0" in "04" and "05" were dropped by tsv export/import
nat_PB_04 = db_nat[["MF"]][db_nat[["SUBSET"]]=="PB" & db_nat[["SUBJECT"]]=="4"]
nat_PB_05 = db_nat[["MF"]][db_nat[["SUBSET"]]=="PB" & db_nat[["SUBJECT"]]=="5"]
nat_PB_11 = db_nat[["MF"]][db_nat[["SUBSET"]]=="PB" & db_nat[["SUBJECT"]]=="11"]

MUT_COL = "MF"

mut_S1 = db[[MUT_COL]][db[["SAMPLE"]]=="S1"]
mut_S1_noMAb = db[[MUT_COL]][db[["SAMPLE"]]=="S1" & !db[["mAb"]]]

mut_S2 = db[[MUT_COL]][db[["SAMPLE"]]=="S2"]
mut_S2_noMAb = db[[MUT_COL]][db[["SAMPLE"]]=="S2" & !db[["mAb"]]]

mut_S3 = db[[MUT_COL]][db[["SAMPLE"]]=="S3"]
mut_S3_noMAb = db[[MUT_COL]][db[["SAMPLE"]]=="S3" & !db[["mAb"]]]

lst = list("05\nNaive"=nat_naive,
           "04\nPB"=nat_PB_04,
           "05\nPB"=nat_PB_05,
           "11\nPB"=nat_PB_11,
           "S1\nPB"=mut_S1,
           "S1\nPB (no mAb)"=mut_S1_noMAb,
           "S2\nPB"=mut_S2,
           "S2\nPB (no mAb)"=mut_S2_noMAb,
           "S3\nPB"=mut_S3,
           "S3\nPB (no mAb)"=mut_S3_noMAb)

#### config ####

# common y-axis
y_max = max(unlist(lst))

cex_count = 0.85
cex_axis_y = 1.25
cex_axis_x = 1.2

col_lst = c("gray50", rep("red", length(lst)-1))

labs_lst = names(lst)

lst_len = unlist(lapply(lst, length))    

#### plot ####

plotName = paste0("mufreq.pdf")
setwd(pathWork)
pdf(plotName, width=12, height=5)
par(mar=c(5,4,0.2,0.2)) # BLTR

beeswarm(lst, method="center", corral="wrap", ylim=c(0, y_max*1.05), 
         xaxt="n", col=col_lst, cex.axis=cex_axis_y, cex=0.5)

for (i in 1:length(lst)) {
    if (length(lst[[i]])>0) {
        boxplot(lst[[i]], outline=F, add=T, at=i, 
                col="transparent", border="gray40", yaxt="n")
    }
}

axis(side=1, at=1:length(lst), line=0, las=2, labels=labs_lst, tick=F, cex.axis=cex_axis_x)
text(x=1:length(lst), y=y_max*1.05, labels=as.character(lst_len), cex=cex_count)

dev.off()

#### special ####

# list of mAbs cross-reactive to human corona viruses
mAbVec = c("98.2B01.H", "98.2E09.H", "98.1E10.H", "05.1C08.H", "88.1A11.H")
# remove ".H"
mAbVec = sapply(mAbVec, function(s){ strsplit(s, "\\.H")[[1]][1] }, USE.NAMES=F)

db_10x = db[!db[["mAb"]], ]

# only 2 out of 5 found to be clonally related to 10x seqs
# TRUE FALSE FALSE FALSE  TRUE
mAbVec %in% db_10x[["CLONAL"]]

# wrt db_10x
bool_cr = db_10x[["CLONAL"]] %in% mAbVec

# 88.1A11 98.2B01 
#    1280      38
table(db_10x[["CLONAL"]][bool_cr], useNA="ifany")

# 10x seqs clonally related cross-reactive mAbs
mut_cr = db_10x[[MUT_COL]][bool_cr]
# rest of 10x seqs 
mut_nonCr = db_10x[[MUT_COL]][!bool_cr]

lst_cr = list("Cross-reactive"=mut_cr, "Rest"=mut_nonCr)

y_max = max(unlist(lst_cr))

cex_count = 0.85
cex_axis_y = 1.25
cex_axis_x = 1.2

col_lst = c("gray50", "gray50")

labs_lst = names(lst_cr)

lst_len = unlist(lapply(lst_cr, length))    

plotName = paste0("mufreq_cr.pdf")
setwd(pathWork)
pdf(plotName, width=5, height=5)
par(mar=c(5,4,0.2,0.2)) # BLTR

beeswarm(lst_cr, method="center", corral="wrap", ylim=c(0, y_max*1.05), 
         xaxt="n", col=col_lst, cex.axis=cex_axis_y, cex=0.5, las=2)

for (i in 1:length(lst_cr)) {
    if (length(lst_cr[[i]])>0) {
        boxplot(lst_cr[[i]], outline=F, add=T, at=i, 
                col="transparent", border="gray40", yaxt="n")
    }
}

axis(side=1, at=1:length(lst_cr), line=0, las=1, labels=labs_lst, tick=F, cex.axis=cex_axis_x)
text(x=1:length(lst_cr), y=y_max*1.05, labels=as.character(lst_len), cex=cex_count)

dev.off()

# 1.473042e-112
wilcox.test(x=mut_cr, y=mut_nonCr)$p.value
