#!/opt/conda/bin/Rscript

# Subsample up to opt$sampleSize number of UIDs with 
# group sizes > opt$LB and < opt$UB

# Output as opt$output with 2 cols: opt$field and SEQCOUNT

suppressPackageStartupMessages(require(optparse))

option_list = list(
    make_option("--input", action="store", default="", type="character", 
                help="Path to input file."),
    make_option("--output", action="store", default="", type="character", 
                help="Path to output file."),
    make_option("--field", action="store", default="", type="character", 
                help="Name of field containing UID."),
    make_option("--LB", action="store", default=20, type="integer", 
                help="Lower bound."),
    make_option("--UB", action="store", default=20000, type="integer", 
                help="Upper bound."),
    make_option("--sampleSize", action="store", default=5000, type="integer", 
                help="Subsample size.")
)
opt = parse_args(OptionParser(option_list=option_list))

# read in ParseHeaders.py table result
db = read.table(opt$input, header=T, sep="\t", stringsAsFactors=F)

# tabulate UIDs
tab = table(db[[opt$field]])

# filter by count
bool = tab > opt$LB & tab < opt$UB 

# subset tab
tab_2 = tab[bool]

# sample from tab_2
n_sample = min(length(tab_2), opt$sampleSize)

set.seed(97458723)
tab_2_sub = tab_2[sample(x=1:length(tab_2), size=n_sample, replace=F)]

# prepare tab_2_sub for export
tab_2_sub_exp = data.frame(matrix(NA, nrow=n_sample, ncol=2))

colnames(tab_2_sub_exp) = c(opt$field, "SEQCOUNT")

tab_2_sub_exp[[opt$field]] = names(tab_2_sub)

tab_2_sub_exp[["SEQCOUNT"]] = as.numeric(tab_2_sub)

write.table(x=tab_2_sub_exp, file=opt$output, sep="\t", 
            quote=F, row.names=F, col.names=T)

cat("Total # sequences selected for EE set:", 
    sum(tab_2_sub_exp[["SEQCOUNT"]]), "\n")
