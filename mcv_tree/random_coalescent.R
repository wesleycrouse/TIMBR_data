library(ape)
trees <- replicate(1000, rcoal(8, LETTERS[1:8]), simplify=F)
write.nexus(trees, file="coalescent_trees.trees")
