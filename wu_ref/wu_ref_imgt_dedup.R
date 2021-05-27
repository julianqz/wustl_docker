# expect a single value from command line input
# version of IMGT reference release
# e.g. "202113-2"
ARGS = commandArgs(trailingOnly=T)
stopifnot(length(ARGS)==1)

#*
PATH_WORK = paste0("~/Dropbox/common/germline_refs/imgt_select/IMGT_vquest_release",
                   ARGS[1])
PATH_SCRIPT = "~/Dropbox/wustl/code/main/shared/bcr_sequence_processing.R"


source(PATH_SCRIPT)

setwd(PATH_WORK)

SPECIES = list.files()
cat("Species:", SPECIES, "\n")

CHAINS = list.files(SPECIES[1])
cat("Chains:", CHAINS, "\n")

for (sp in SPECIES) {
    for (ch in CHAINS) {
        
        all_files = list.files(path=paste(".", sp, ch, sep="/"),
                               pattern="fasta")
        
        for (f in all_files) {
            
            # read in multi-line fasta file
            path_f = paste0(paste(".", sp, ch, sep="/"), "/", f)
            vec = remove_duplicate_fasta(path_f)
            
            # parse fasta filename
            # get the part in front of ".fasta"
            filename_stem = strsplit(f, split="\\.")[[1]][1]
            
            # filename for exporting
            filename_export = paste0(filename_stem, "_no_dup.fasta")
            
            # path for export file
            path_export = paste0(paste(".", sp, ch, sep="/"), "/",
                                 filename_export)
            # export
            # add_header_symbol=F because `read_multiline_fasta` keeps `>` in headers
            cat("Exported:", filename_export, "\n")
            export_fasta(sequences=vec, headers=names(vec),
                         add_header_symbol=FALSE, filename=path_export)
            
        }
        
    }
}
