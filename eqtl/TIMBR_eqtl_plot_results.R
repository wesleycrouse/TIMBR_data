load("eqtl_full_results.RData")

#calculate expected number of alleles and sort
results.table <- t(sapply(output, function(x){c(x$probe, sum(x$p.K.given.y*as.numeric(names(x$p.K.given.y))))}))
results.table <- as.data.frame(results.table, stringsAsFactors=F)
colnames(results.table) <- c("Probe", "Exp.K")
results.table$Exp.K <- as.numeric(results.table$Exp.K)
results.table <- results.table[order(-results.table$Exp.K),]

#merge with probe names
eqtl.list <- read.csv("data/eQTL_list.csv", stringsAsFactors=F)
eqtl.list <- eqtl.list[eqtl.list$category %in% c("cis", "cis(qtl_in_gene)", "near-cis") & eqtl.list$merge_SNP==1,c(1,2,4,5)]
eqtl.list[,1] <- sapply(eqtl.list[,1], function(x){paste0("X", substring(x,2))})
eqtl.list <- as.data.frame(t(sapply(results.table$Probe, function(x){eqtl.list[which(x==eqtl.list),]})))

results.table <- cbind(results.table, eqtl.list[,c(2,3,4)])
results.table <- results.table[,c(1,3:5,2)]
colnames(results.table)[2:4] <- c("Gene", "Chr", "Pos")
rownames(results.table) <- NULL

#table of ~1% most multiallelic genes
head(results.table, floor(nrow(results.table)*0.01))

library(xtable)
table.data <- head(results.table, 20)
table.data <- table.data[,c(1,2,3,5)]
xtable(table.data, digits=c(0,0,0,0,4))

#distribution of number of alleles over all genes
full.distribution <- unlist(sapply(output, function(x){x$p.K.given.y}))
full.distribution <- tapply(full.distribution, names(full.distribution), sum)/length(output)

prior.alleles <- exp(sapply(1:8, TIMBR:::ln.K.prior.crp.marginalized, J=8, a=1, b=2.333415))

png(filename = "fig_eqtl_alleles.png", width = 600, height = 600)

df.plot <- barplot(full.distribution, xlab="Alleles", ylab="Posterior Probability", main="Distribution of Number of Alleles", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, ylim=c(0,1))
lines(df.plot, prior.alleles, lty=5)
points(df.plot, prior.alleles, pch=15)

dev.off()

####################
#Glo1 results
library(TIMBR)

load("results/results_crp_X2667352.RData")

results.full <- TIMBR(results.crp$y, results.crp$prior.D, list(model.type="fixed", M.IDs=paste(0:7, collapse=",")), samples=100000)
save(results.full, file="results/results_crp_X2667352_full.RData")
#load("results/results_crp_X2667352_full.RData")

#Bayes factor
results.crp$ln.BF - results.full$ln.BF

#top allelic series
head(results.crp$p.M.given.y, 10)

#distribution of number of alleles for Glo1
png(filename = "fig_glo1_alleles.png", width = 600, height = 600)

posterior.alleles <- c(rep(0,3), results.crp$p.K.given.y)
names(posterior.alleles)[1:3] <- 1:3

df.plot <- barplot(posterior.alleles, xlab="Alleles", ylab="Posterior Probability", main="Distribution of Number of Alleles", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, ylim=c(0,1))
lines(df.plot, prior.alleles, lty=5)
points(df.plot, prior.alleles, pch=15)

dev.off()

#distribution of haplotype effects for Glo1

png(filename = "fig_glo1_haplotypes.png", width = 480, height = 480)

TIMBR.plot.haplotypes(results.crp, TIMBR.output.bkgrd=results.full, x.lab="Glo1", x.lim=c(-2.5, 2.5))
title(main="CRP vs Full Model")

dev.off()

#plot raw phenotype

png(filename = "fig_glo1_phenotype.png", width = 480, height = 480)

D <- apply(results.crp$prior.D$P, 1, which.max)
y <- results.crp$y[D<=8]
D <- -D[D<=8]

colors <- rev(c("#9000E0","#F00000","#00A000","#00A0F0","#1010F0","#F08080","#808080","#F0F000"))

plot(y, jitter(D, 0.5), bg=colors[-D], pch=21, xlab="Glo1", ylab="Haplotype", cex=1.3, yaxt='n', main="Observed Phenotype", xlim=c(-2.5, 2.5))
axis(side=2, at=-1:-8, labels=LETTERS[1:8], tick=F, las=1)

dev.off()