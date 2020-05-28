eqtl.list <- read.csv("data/eQTL_list.csv", stringsAsFactors=F)
eqtl.list <- eqtl.list[eqtl.list$category %in% c("cis", "cis(qtl_in_gene)", "near-cis") & eqtl.list$merge_SNP==1, c("probe","peak_locus")]
eqtl.list[,1] <- sapply(eqtl.list[,1], function(x){paste0("X", substring(x,2))})

write.table(eqtl.list, file="eqtl_list.txt", quote=F, row.names=F, col.names=F)
