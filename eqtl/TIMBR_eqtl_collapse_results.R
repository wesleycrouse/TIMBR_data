file.list <- list.files("results/")

output <- list()

for (i in file.list){
	probe <- unlist(strsplit(unlist(strsplit(i, "_"))[3], "[.]"))[1]
	print(probe)
	load(paste0("results/", i))
	output[[length(output)+1]] <- list(probe=probe, p.K.given.y=results.crp$p.K.given.y)
}

save(output, file="eqtl_full_results.RData")
