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
prior.M <- list(model.type="fixed", M.IDs=paste(0:7, collapse=","))
results.full <- TIMBR(y, prior.D, prior.M, samples=100000)
save(results.full, file="results_full_mcv.RData")
#load("results_full_mcv.RData")

#call TIMBR with CRP model
prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=1, prior.alpha.rate=2.333415)
results.crp <- TIMBR(y, prior.D, prior.M, samples=100000)
save(results.crp, file="results_crp_mcv.RData")
#load("results_crp_mcv.RData")

#call TIMBR with tree model
load("../mcv_tree/mcv_trees.RData")
trees.prior <- lapply(trees, ewenss.calc, prior.alpha=list(type="gamma", shape=1, rate=2.333415))
save(trees.prior, file="mcv_trees_prior.RData")
#load("mcv_trees_prior.RData")

trees.prior <- unlist(lapply(trees.prior, function(x){ln.probs <- x$ln.probs; names(ln.probs) <- x$M.IDs; ln.probs}))
trees.prior <- tapply(trees.prior, names(trees.prior), matrixStats::logSumExp)
prior.M <- list(model.type="list", M.IDs=names(trees.prior), ln.probs=as.numeric(unname(trees.prior)-log(1000)), hash.names=T)

results.tree <- TIMBR(y, prior.D, prior.M, samples=100000)
save(results.tree, file="results_tree_mcv.RData")
#load("results_tree_mcv.RData")

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
#generate plots
#prior probabilities jitter plot
prior.M.crp <- ewenss.calc(8, list(type="gamma", shape=1, rate=2.333415))

png(filename = "fig_mcv_prior_jitter.png", width = 960, height = 480)

par(mfrow = c(1, 2))

plot(jitter(sapply(prior.M.crp$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1})), prior.M.crp$ln.probs,
     ylab="ln(Prior Probability)", xlab="Number of Alleles", ylim=c(-25, 0), xlim=c(1,8), main="CRP")

plot(jitter(sapply(prior.M$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1})), prior.M$ln.probs,
     ylab="ln(Prior Probability)", xlab="Number of Alleles", ylim=c(-25, 0), xlim=c(1,8), main="Tree")

dev.off()

#distribution of posterior number of alleles for CRP and tree

png(filename = "fig_mcv_alleles.png", width = 600, height = 600)

allele.table <- table(c(results.crp$post.K, results.tree$post.K), c(rep(0, length(results.crp$post.K)), rep(1, length(results.tree$post.K))))
allele.table <- rbind(0, allele.table)/results.tree$samples
row.names(allele.table)[1] <- 1
df.plot <- barplot(t(allele.table), beside=T, xlab="Alleles", ylab="Posterior Probability", main="Distribution of Number of Alleles", 
                   cex.lab=1.5, cex.axis=1.2, cex.main=1.5, legend=c("CRP", "Tree"), ylim=c(0, 1))


prior.alleles <- tapply(exp(prior.M$ln.probs), sapply(prior.M$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1}), sum)
prior.alleles.crp <- tapply(exp(prior.M.crp$ln.probs), sapply(prior.M.crp$M.IDs, function(x){max(as.numeric(unlist(strsplit(x, ","))))+1}), sum)

lines(df.plot[1,]+0.5, prior.alleles.crp, lty=5)
points(df.plot[1,]+0.5, prior.alleles.crp, pch=15)

lines(df.plot[1,]+0.5, prior.alleles, lty=3)
points(df.plot[1,]+0.5, prior.alleles, pch=16)

dev.off()

#distribution of haplotype effects for CRP and tree

png(filename = "fig_mcv_haplotypes.png", width = 960, height = 480)

par(mfrow = c(1, 2))

TIMBR.plot.haplotypes(results.crp, TIMBR.output.bkgrd=results.full, x.lim=c(35, 75), x.lab="Mean Cell Volume (fL)")
title(main="CRP vs Full Model")

TIMBR.plot.haplotypes(results.tree, TIMBR.output.bkgrd=results.crp, x.lim=c(35, 75), x.lab="Mean Cell Volume (fL)", colors=scales::alpha(rep("#E6E6E6", 8), 0.4), colors.bkgrd=scales::alpha(rep("#4D4D4D", 8), 0.7))
title(main="Tree vs CRP Model")

dev.off()

#plot raw phenotype

png(filename = "fig_mcv_phenotype.png", width = 480, height = 480)

D <- apply(results.crp$prior.D$P, 1, which.max)
y <- results.crp$y[D<=8]
D <- -D[D<=8]

colors <- rev(c("#9000E0","#F00000","#00A000","#00A0F0","#1010F0","#F08080","#808080","#F0F000"))

plot(y, jitter(D, 0.5), bg=colors[-D], pch=21, xlab="Mean Cell Volume (fL)", ylab="Haplotype", cex=1.3, yaxt='n', main="Observed Phenotype", xlim=c(35, 75))
axis(side=2, at=-1:-8, labels=LETTERS[1:8], tick=F, las=1)

dev.off()
