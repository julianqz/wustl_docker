# Julian Q. Zhou
# https://github.com/julianqz

# De-duplicate IMGT Ig allele references

# expect a single value from command line input
# version of IMGT reference release

# e.g. "202113-2"
ARGS = commandArgs(trailingOnly=T)
stopifnot(length(ARGS)==1)

#*
PATH_WORK_1 = "~/Dropbox/common/germline_refs/imgt_select/"
#PATH_WORK_1 = "~/Dropbox/common/germline_refs/C57BL6/"
#PATH_WORK_1 = "~/Dropbox/common/germline_refs/H2L2/"

PATH_WORK_2 = paste0(PATH_WORK_1, "IMGT_vquest_release", ARGS[1])
#PATH_WORK_2 = paste0(PATH_WORK_1, "IMGT_genedb_", ARGS[1]) #* H2L2

PATH_SCRIPT = "~/Dropbox/wustl/code/c2b2/bcr_sequence_processing.R"


source(PATH_SCRIPT)

setwd(PATH_WORK_2)

SPECIES = list.files()
cat("Species:", SPECIES, "\n")

CHAINS = list.files(SPECIES[1])
cat("Chains:", CHAINS, "\n")

for (sp in SPECIES) {
    for (ch in CHAINS) {

        #### individual files ####
        
        all_files = list.files(path=paste(".", sp, ch, sep="/"),
                               pattern="fasta")
        
        for (f in all_files) {
            
            # read in multi-line fasta file
            path_f = paste0(paste(".", sp, ch, sep="/"), "/", f)
            vec = remove_duplicate_fasta(path_f)
            
            #*
            RUN_MOD_IMGT_ACCESSION=F
            
            if (RUN_MOD_IMGT_ACCESSION) {
                
                #* IMGT-specific
                # headers: 
                # >AC003996|TRAV3-4*01|Mus musculus_129/SvJ|F|V-REGION|30724..31003|280 nt|1| | | | |280+48=328| | |
                # >AC003996|TRAV13-4/DV7*01|Mus musculus_129/SvJ|F|V-REGION|4043..4313|271 nt|1| | | | |271+57=328| | |
                # Notice how the accession number after > can be the same for entries with unique nt seqs
                # De-dup accession numbers
                
                # accession number
                vec_imgt_header_acc = sapply(names(vec), 
                                             function(s){strsplit(s, "\\|")[[1]][1]}, 
                                             USE.NAMES=F)
                # rest of header
                vec_imgt_header_rest = sapply(names(vec), 
                                              function(s){strsplit(s, "^>\\w+\\|")[[1]][2]},
                                              USE.NAMES=F)
                
                stopifnot(all.equal( paste0(vec_imgt_header_acc, "|", vec_imgt_header_rest),
                                     names(vec) ))
                
                if (any(duplicated(vec_imgt_header_acc))) {
                    tab_acc = table(vec_imgt_header_acc)
                    # unique accession numbers that have duplicates
                    vec_acc_dup_uniq = names(tab_acc[tab_acc>1])
                    
                    # for each unique accession number, de-dup by adding suffix _[123..]
                    for (cur_uniq in vec_acc_dup_uniq) {
                        # wrt vec_imgt_header_acc
                        cur_uniq_idx = which(vec_imgt_header_acc==cur_uniq)
                        # add suffix
                        cur_dedup_acc = paste0(cur_uniq, "_", 
                                               as.character(1:length(cur_uniq_idx)))
                        
                        cat("- de-dup IMGT accession number in header from",
                            cur_uniq, "into:\n", " ", cur_dedup_acc, "\n")
                        vec_imgt_header_acc[cur_uniq_idx] = cur_dedup_acc
                    }
                    
                    # reset header
                    vec_imgt_header_dedup = paste0(vec_imgt_header_acc, "|", vec_imgt_header_rest)
                    stopifnot(!any(duplicated(vec_imgt_header_dedup)))
                    
                    names(vec) = vec_imgt_header_dedup
                } 
                
            }

            
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
        
        #### by segment ####
        
        # as of 2025-05-23, only done for TR
        
        if (ch=="TR") {
            
            vec_no_dup_files = list.files(path=paste(".", sp, ch, sep="/"),
                                          pattern="_no_dup.fasta")
            
            for (seg in c("V", "D", "J", "C")) {
                
                cat("segment:", seg, "\n")
                vec_no_dup_seg_files = vec_no_dup_files[grepl(pattern=paste0("", seg, "_no_dup.fasta"), 
                                                              x=vec_no_dup_files, fixed=T)]
                
                lst_no_dup_seg = vector(mode="list", length=length(vec_no_dup_seg_files))
                
                for (i_file in 1:length(vec_no_dup_seg_files)) {
                    
                    i_f = vec_no_dup_seg_files[i_file]
                    cat(i_f, "\n")
                    
                    i_file_name = paste0(paste(".", sp, ch, sep="/"), "/", i_f)
                    lst_no_dup_seg[[i_file]] = read_multiline_fasta(i_file_name)
                }
                
                vec_no_dup_seg_all = do.call(c, lst_no_dup_seg)
                
                vec_no_dup_seg_all_dedup = remove_duplicate_fasta(vec_fasta=vec_no_dup_seg_all)
                
                i_fn_export = paste0("concat_no_dup_", seg, ".fasta")
                i_path_export = paste0(paste(".", sp, ch, sep="/"), "/", i_fn_export)
                
                cat("Exported:", i_path_export, "\n")
                export_fasta(sequences=vec_no_dup_seg_all_dedup, 
                             headers=names(vec_no_dup_seg_all_dedup),
                             add_header_symbol=F,
                             filename=i_path_export)
            }
            
        }
   
    }
}
