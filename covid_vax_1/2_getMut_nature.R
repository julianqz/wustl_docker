pathWork = "~/Dropbox (recherche)/wustl/projects/covid_vax_1/JZ/mutation/"

RUN=F
if (RUN) {
    
    # Manually run bcr_8_muFreq_plot_10x.R
    # Then save lst_plot[["Naive\n(PBMC)"]]
    # P05 total 10x naive
    # Extended Data Fig 3h, naive from PBMC (Turner et al 2020)
    
    nat_naive = lst_plot[["Naive\n(PBMC)"]]
    length(nat_naive) # 2553
    
    fn = "nat_naive.RData"
    setwd(pathWork)
    save(nat_naive, file=fn)
    
    # Manually run bar_8_muFreq_plot_bulk.R
    # Then save lst_plot[[1]]; lst_plot[[4]]; lst_plot[[7]]
    # P04, P05, P11's "nonShared_PB"
    # bool_mAb_clone & db[["shared_GC_PB"]]=="N" & db[["CELL_TYPE"]]=="PB"
    # Fig 3c (Turner et al 2020)
    
    nat_PB_04 = lst_plot[[1]]
    nat_PB_05 = lst_plot[[4]]
    nat_PB_11 = lst_plot[[7]]
    
    length(nat_PB_04) # 149
    length(nat_PB_05) # 1033
    length(nat_PB_11) # 56
    
    fn = "nat_PB.RData"
    setwd(pathWork)
    save(nat_PB_04, nat_PB_05, nat_PB_11, file=fn)
    
    ### export
    
    db_naive = data.frame(cbind("SUBJECT"="05", "SUBSET"="Naive", "MF"=nat_naive))
    db_PB_04 = data.frame(cbind("SUBJECT"="04", "SUBSET"="PB", "MF"=nat_PB_04))
    db_PB_05 = data.frame(cbind("SUBJECT"="05", "SUBSET"="PB", "MF"=nat_PB_05))
    db_PB_11 = data.frame(cbind("SUBJECT"="11", "SUBSET"="PB", "MF"=nat_PB_11))
    db_nat = rbind(db_naive, db_PB_04, db_PB_05, db_PB_11)
    
    # sanity check
    #    Naive   PB
    # 04     0  149
    # 05  2553 1033
    # 11     0   56
    table(db_nat[["SUBJECT"]], db_nat[["SUBSET"]])
    
    fn_nat = "mufreq_nature.tsv"
    setwd(pathWork)
    write.table(x=db_nat, file=fn_nat, quote=F, sep="\t",
                row.names=F, col.names=T)
}

