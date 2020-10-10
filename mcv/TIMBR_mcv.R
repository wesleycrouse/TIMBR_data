library(TIMBR)

#load phenotype data
phenotype.data <- read.table("data/mcv_pheno.txt", header=T, stringsAsFactors=F)

#drop subjects with any missing phenotype data
phenotype.data <- phenotype.data[!is.na(phenotype.data[,2]),]

#sort phenotype data by subject name
phenotype.data <- phenotype.data[order(phenotype.data[,1]),]

#format dependent variable
y <- phenotype.data[,2]
names(y) <- phenotype.data[,1]

#load diplotype data
locus <- "M7.961"

chr <- substring(unlist(strsplit(locus, split="[.]"))[1], 2)
load(paste0("../precc_cache/full/chr", chr, "/@", locus, ".RData"))
load(paste0("../precc_cache/full/chr", chr, "/subjects.RData"))

P <- get(locus)
rownames(P) <- subjects

do.call("rm", list(eval(locus)))
rm(subjects)

#sort diplotype probabilities by subject name
P <- P[order(rownames(P)),]

#drop subjects without phenotype data
P <- P[rownames(P) %in% phenotype.data[,1],]

#specify prior.D
prior.D <- list(P=P, A=additive.design(8, "happy"), fixed.diplo=F)

#call TIMBR with full model
#prior.M <- list(model.type="fixed", M.IDs=paste(0:7, collapse=","))
#results.full <- TIMBR(y, prior.D, prior.M, samples=100000)
#save(results.full, file="results_full_mcv.RData")
load("results_full_mcv.RData")

#call TIMBR with CRP model
#prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=1, prior.alpha.rate=2.333415)
#results.crp <- TIMBR(y, prior.D, prior.M, samples=100000)
#save(results.crp, file="results_crp_mcv.RData")
load("results_crp_mcv.RData")

#call TIMBR with tree model
load("../mcv_tree/mcv_trees.RData")
#trees.prior <- lapply(trees, ewenss.calc, prior.alpha=list(type="gamma", shape=1, rate=2.333415))
#save(trees.prior, file="mcv_trees_prior.RData")
load("mcv_trees_prior.RData")

trees.prior <- unlist(lapply(trees.prior, function(x){ln.probs <- x$ln.probs; names(ln.probs) <- x$M.IDs; ln.probs}))
trees.prior <- tapply(trees.prior, names(trees.prior), matrixStats::logSumExp)
prior.M <- list(model.type="list", M.IDs=names(trees.prior), ln.probs=as.numeric(unname(trees.prior)-log(1000)), hash.names=T)

#results.tree <- TIMBR(y, prior.D, prior.M, samples=100000)
#save(results.tree, file="results_tree_mcv.RData")
load("results_tree_mcv.RData")

####################
#report statistics
#calculate median, HPD and length of interval for each approach
library(coda)
library(xtable)

colnames(results.full$post.hap.effects) <- LETTERS[1:8]
colnames(results.crp$post.hap.effects) <- LETTERS[1:8]
colnames(results.tree$post.hap.effects) <- LETTERS[1:8]

HPD.full <- HPDinterval(as.mcmc(results.full$post.hap.effects))
HPD.crp <- HPDinterval(as.mcmc(results.crp$post.hap.effects))
HPD.tree <- HPDinterval(as.mcmc(results.tree$post.hap.effects))

#median
cbind(apply(results.full$post.hap.effects, 2, median), 
      apply(results.crp$post.hap.effects, 2, median), 
      apply(results.tree$post.hap.effects, 2, median))

#HPD
HPD.full
HPD.crp
HPD.tree

#HPD length
xtable(cbind(apply(HPD.full, 1, function(x){x[2]-x[1]}), 
             apply(HPD.crp, 1, function(x){x[2]-x[1]}),
             apply(HPD.tree, 1, function(x){x[2]-x[1]})), digits=2)

cbind(apply(HPD.full, 1, function(x){x[2]-x[1]}), 
      apply(HPD.crp, 1, function(x){x[2]-x[1]}),
      apply(HPD.tree, 1, function(x){x[2]-x[1]}))

#top allelic series
head(results.crp$p.M.given.y, 10)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")

table.data <- head(results.crp$p.M.given.y, 10)
table.data <- cbind(sapply(names(table.data), function(x){max(TIMBR:::m.from.M.ID(x))}), table.data)
xtable(table.data, digits=c(0,0,4))

table.data <- head(results.tree$p.M.given.y, 10)
table.data <- cbind(sapply(names(table.data), function(x){max(TIMBR:::m.from.M.ID(x))}), table.data)
xtable(table.data, digits=c(0,0,4))

table.data <- exp(prior.M$ln.probs)
names(table.data) <- prior.M$M.IDs
table.data <- -sort(-table.data)
table.data <- cbind(sapply(names(table.data), function(x){max(TIMBR:::m.from.M.ID(x))}), table.data)
xtable(table.data[1:10,], digits=c(0,0,4))

head(results.tree$p.M.given.y, 10)

results.crp$p.M.given.y[results.crp$p.M.given.y>0.01]
results.tree$p.M.given.y[results.tree$p.M.given.y>0.01]

sum(results.crp$p.M.given.y[results.crp$p.M.given.y>0.01])
sum(results.tree$p.M.given.y[results.tree$p.M.given.y>0.01])
sum(results.crp$p.M.given.y[results.crp$p.M.given.y>0.01])-results.crp$p.M.given.y[1]

#expectation of number of alleles
mean(results.crp$post.K)
mean(results.tree$post.K)

#bayes factor
results.crp$ln.BF - results.full$ln.BF
results.tree$ln.BF - results.crp$ln.BF
results.tree$ln.BF - results.full$ln.BF

#number of topologies
#vs 10395 total topologies: https://en.wikipedia.org/wiki/Unrooted_binary_tree
trees.unique <- ape::unique.multiPhylo(trees)
length(trees.unique)
table(attr(trees.unique, "old.index"))

#number of prior allelic series
length(prior.M$M.IDs)

####################
#function for legend
source("../plot_functions/plot_functions.R")

####################
#generate plots
#distribution of posterior number of alleles for CRP and tree
prior.M.crp <- ewenss.calc(8, list(type="gamma", shape=1, rate=2.333415))

allele.table <- table(c(results.crp$post.K, results.tree$post.K), c(rep(0, length(results.crp$post.K)), rep(1, length(results.tree$post.K))))
allele.table <- rbind(0, allele.table)/results.tree$samples
row.names(allele.table)[1] <- 1

colors <- c(scales::alpha("orange", 0.8),
            scales::alpha("blue", 0.5))

pdf(file = "fig_mcv_alleles.pdf", width = 4, height = 4)
par(mai=c(0.75, 0.8, 0.5, 0.2))

barplot(matrix(NA, 2, 8), 
        beside=T, 
        xlab="", 
        ylab="Probability", 
        main="Distribution of Alleles", 
        font.main=1,
        las=1,
        cex.lab=1.1, 
        cex.axis=1, 
        cex.main=1.2, 
        ylim=c(0, 1))

title(xlab="Number of Alleles", line=2.5, cex.lab=1.1)

box()

for (h in (1:5)/5){
   abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
}

df.plot <- barplot(t(allele.table), beside=T, 
                   yaxt='n',
                   col=colors, add=T)

legend_homebrew("topright", inset=.02, legend=c("CRP Posterior", "Tree Posterior", "CRP Prior", "Tree Prior"),
       fill=c(colors, NA, NA),
       border=c("black", "black", "white", "white"),
       bg="white",
       cex=0.8,
       pch=c(NA,NA,22,22),
       pt.bg=c(NA, NA, "orange", "blue"),
       lwd=1,
       lty=c(NA,NA,5,5),
       col=c(NA, NA, "orange", "blue"))

prior.alleles <- tapply(exp(prior.M$ln.probs), sapply(prior.M$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1}), sum)
prior.alleles.crp <- tapply(exp(prior.M.crp$ln.probs), sapply(prior.M.crp$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1}), sum)

lines(df.plot[1,]+0.5, prior.alleles.crp, lty=5, lwd=1, col="orange")
points(df.plot[1,]+0.5, prior.alleles.crp, pch=22, cex=1, bg="orange")

lines(df.plot[1,]+0.5, prior.alleles, lty=5, lwd=1, col="blue")
points(df.plot[1,]+0.5, prior.alleles, pch=22, cex=1, bg="blue")

dev.off()

####################
#combined phenotype and haplotype plot

pdf(file = "fig_mcv_phenohap.pdf", width = 8, height = 4)

par(mfrow = c(1, 2))
par(mai=c(0.8, 0.8, 0.5, 0.2))

D <- apply(results.crp$prior.D$P, 1, which.max)
y <- results.crp$y[D<=8]
D <- -D[D<=8]

plot(y, jitter(D, 0.3),
     bg=scales::alpha(rep("#4D4D4D", 8), 0.7),
     pch=21, 
     xlab="", 
     ylab="Haplotype",
     main="Observed Values", 
     font.main=1,
     cex=1.2, 
     yaxt='n', 
     xlim=c(35, 75),
     cex.lab=1, 
     cex.axis=1, 
     cex.main=1.1)

axis(side=2, at=-1:-8, labels=LETTERS[1:8], tick=F, las=1)
title(xlab="MCV (fL)", line=2.5)

put.fig.letter("A", "topleft", font=2, cex=1.2, offset=c(0.005,-0.005))

TIMBR.plot.haplotypes(results.crp, TIMBR.output.bkgrd=results.full, 
                      x.lim=c(35, 75),
                      x.lab="",
                      cex.lab=1, 
                      cex.axis=1, 
                      cex.main=1.1, 
                      font.main=1,
                      colors=scales::alpha(rep("orange", 8), 0.8),
                      colors.bkgrd=scales::alpha(rep("#4D4D4D", 8), 0.8))

title(xlab="MCV (fL)", line=2.5)
title(main="Modeled Effects")

box()

legend("topleft", 
       inset=.02, 
       c("CRP", "Full"), 
       fill=c(scales::alpha("orange", 0.8), 
              scales::alpha("#4D4D4D", 0.8)),
       cex=0.8)

put.fig.letter("B", "topleft", font=2, cex=1.2, offset=c(0.003,0))

dev.off()

####################
#prior probabilities jitter plot
pdf(file = "fig_mcv_prior.pdf", width = 8, height = 8)

par(mfrow = c(2, 2))
par(mai=c(0.7, 0.8, 0.7, 0.25))

plot(c(), c(), 
     ylim=c(1,8), 
     xlim=c(1,8),
     cex.lab=1.2, 
     cex.axis=1.1, 
     cex.main=1.3, 
     font.main = 1,
     xaxt="n", yaxt="n", ann=F)

title(main="Coalescent", cex.main=1.3, font.main = 1, line=1)
title(xlab="Time", cex.lab=1.2, font.main = 1, line=1)
title(ylab="Haplotype", cex.lab=1.2, font.main = 1, line=3)
axis(2, at=8:1, labels=c("A","D","G","F","C","B","H","E"), las=1, tick=T)
axis(4, at=8:1, labels=F, las=1, tick=T, tcl=0.5)

put.fig.letter("A", "topleft", font=2, cex=1.5, offset=c(0.005, -0.06))

plot(c(), c(), 
     ylim=c(1,8), 
     xlim=c(1,8),
     cex.lab=1.2, 
     cex.axis=1.1, 
     cex.main=1.3, 
     font.main = 1,
     xaxt="n", yaxt="n", ann=F)

title(main="MCV Locus", cex.main=1.3, font.main = 1, line=1)
title(xlab="Time", cex.lab=1.2, font.main = 1, line=1)
title(ylab="Haplotype", cex.lab=1.2, font.main = 1, line=3)
axis(2, at=8:1, labels=c("B","E","D","G","F","C","H","A"), las=1, tick=T)
axis(4, at=8:1, labels=F, las=1, tick=T, tcl=0.5)

put.fig.letter("B", "topleft", font=2, cex=1.5, offset=c(0.005, -0.06))

par(mai=c(0.8, 0.8, 0.55, 0.25))

plot(c(), c(), 
     ylab="Prior Probability (ln)", 
     xlab="Number of Alleles", 
     ylim=c(-25, 0), 
     xlim=c(1,8), 
     main="CRP",
     cex.lab=1.2, 
     cex.axis=1.1, 
     cex.main=1.3, 
     font.main = 1, 
     las=1)

for (h in -5*(0:5)){
   abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.8)
}

points(jitter(sapply(prior.M.crp$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1})), 
       prior.M.crp$ln.probs,
       lwd=1,
       pch=21,
       bg=scales::alpha("orange", 0.8),
       cex=1.3)

put.fig.letter("C", "topleft", font=2, cex=1.5, offset=c(0.005, -0.005))

plot(c(), c(), 
     ylab="Prior Probability (ln)", 
     xlab="Number of Alleles", 
     ylim=c(-25, 0), 
     xlim=c(1,8), 
     main="Tree",
     cex.lab=1.2, 
     cex.axis=1.1, 
     cex.main=1.3, 
     font.main = 1, 
     las=1)

for (h in -5*(0:5)){
   abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.8)
}

points(jitter(sapply(prior.M$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1})), 
       prior.M$ln.probs,
       lwd=1,
       pch=21,
       bg=scales::alpha("blue", 0.8),
       cex=1.1)

put.fig.letter("D", "topleft", font=2, cex=1.5, offset=c(0.005, -0.005))

dev.off()
