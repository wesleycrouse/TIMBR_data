library(ape)

png(filename = "tree1.png", width = 480, height = 480)

set.seed(20)
tree1 <- rcoal(4, tip.label=LETTERS[1:4])
plot.phylo(tree1, type="cladogram", edge.width=3, font=2, label.offset=0.02)

dev.off()

png(filename = "tree2.png", width = 480, height = 480)

set.seed(14)
tree2 <- rcoal(4, tip.label=LETTERS[1:4])
plot.phylo(tree2, type="cladogram", edge.width=3, font=2, label.offset=0.02)

dev.off()