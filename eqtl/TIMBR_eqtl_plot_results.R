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

####################
source("../plot_functions/plot_functions.R")

####################
#distribution of number of alleles over all genes

full.distribution <- unlist(sapply(output, function(x){x$p.K.given.y}))
full.distribution <- tapply(full.distribution, names(full.distribution), sum)/length(output)

prior.alleles <- exp(sapply(1:8, TIMBR:::ln.K.prior.crp.marginalized, J=8, a=1, b=2.333415))

pdf(file = "fig_eqtl_alleles.pdf", width = 4, height = 4)

par(mai=c(0.75, 0.8, 0.5, 0.2))

df.plot <- barplot(rep(NA,8), 
                   xlab="", 
                   ylab="Probability",
                   main="Average Distribution of Alleles", 
                   font.main=1,
                   las=1,
                   cex.lab=1.1, 
                   cex.axis=1, 
                   cex.main=1.2, 
                   ylim=c(0,1))

title(xlab="Number of Alleles", line=2.5, cex.lab=1.1)

box()

for (h in (1:5)/5){
  abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
}

df.plot <- barplot(full.distribution, add=T,
                   yaxt='n',
                   col=scales::alpha("orange", 0.8))

lines(df.plot, prior.alleles, lty=5, lwd=1, col="orange")
points(df.plot, prior.alleles, pch=22, cex=1, bg="orange")

legend_homebrew("topright", inset=.02, legend=c("CRP Posterior", "CRP Prior"),
                fill=c(scales::alpha("orange", 0.8), NA),
                border=c("black", "white"),
                bg="white",
                cex=0.8,
                pch=c(NA,22),
                pt.bg=c(NA, "orange"),
                lwd=1,
                lty=c(NA,5),
                col=c(NA, "orange"))

dev.off()

####################
#Glo1 results
library(TIMBR)

load("results/results_crp_X2667352.RData")

#results.full <- TIMBR(results.crp$y, results.crp$prior.D, list(model.type="fixed", M.IDs=paste(0:7, collapse=",")), samples=100000)
#save(results.full, file="results/results_crp_X2667352_full.RData")
load("results/results_crp_X2667352_full.RData")

#Bayes factor
results.crp$ln.BF - results.full$ln.BF

#top allelic series
head(results.crp$p.M.given.y, 10)

####################
#raw phenotype and distribution of haplotype effects for Glo1

pdf(file = "fig_glo1.pdf", width = 8, height = 4)

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
     xlim=c(-2.5, 2.5),
     cex.lab=1, 
     cex.axis=1, 
     cex.main=1.1)

axis(side=2, at=-1:-8, labels=LETTERS[1:8], tick=F, las=1)
title(xlab="Glo1", line=2.5)

put.fig.letter("A", "topleft", font=2, cex=1.2, offset=c(0.005,-0.005))

posterior.alleles <- c(rep(0,3), results.crp$p.K.given.y)
names(posterior.alleles)[1:3] <- 1:3

df.plot <- barplot(rep(NA,8), 
                   xlab="", 
                   ylab="Probability",
                   main="Distribution of Alleles", 
                   font.main=1,
                   las=1,
                   cex.lab=1,
                   cex.axis=1, 
                   cex.main=1.1, 
                   ylim=c(0,1))

title(xlab="Alleles", line=2.5)

box()

for (h in (1:5)/5){
  abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
}

df.plot <- barplot(posterior.alleles, 
                   add=T,
                   yaxt='n',
                   col=scales::alpha("orange", 0.8))

lines(df.plot, prior.alleles, lty=5, lwd=1, col="orange")
points(df.plot, prior.alleles, pch=22, cex=1, bg="orange")

legend_homebrew("topright", inset=.02, legend=c("CRP Posterior", "CRP Prior"),
                fill=c(scales::alpha("orange", 0.8), NA),
                border=c("black", "white"),
                bg="white",
                cex=0.8,
                pch=c(NA,22),
                pt.bg=c(NA, "orange"),
                lwd=1,
                lty=c(NA,5),
                col=c(NA, "orange"))

put.fig.letter("B", "topleft", font=2, cex=1.2, offset=c(0.003,0))

dev.off()
