library(TIMBR)
library(ape)

burn.in <- 101

#load trees from BEAST and discard burn-in
trees <- read.nexus("chr7_3776.trees")
trees <- trees[-c(1:burn.in)]

#load log file and extract population sizes
log.file <- read.table("chr7_3776.log", skip=3, header = TRUE)
pop.sizes <- log.file$constant.popSize[-c(1:burn.in)]
rm(log.file, burn.in)

#convert edge lengths to coalescent units
trees <- lapply(1:length(trees), function(x){trees[[x]]$edge.length <- trees[[x]]$edge.length/pop.sizes[x]; trees[[x]]})
rm(pop.sizes)

save(trees, file="mcv_trees.RData")
write.nexus(trees, file="mcv_trees.trees")
