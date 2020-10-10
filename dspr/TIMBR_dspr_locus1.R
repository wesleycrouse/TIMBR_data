phenotype.var <- "CG10245"
cache.file <- "P_locus1.txt"

####################
#devtools::install_github("wesleycrouse/TIMBR")
library(TIMBR)

#load phenotype data
phenotype.data <- read.table(file="data/FemaleHeadExpression_trimmed.txt", header=T, sep=" ")
phenotype.data <- cbind(phenotype.data, NA)
colnames(phenotype.data)[ncol(phenotype.data)] <- "SUBJECT.NAME"

for (i in 1:nrow(phenotype.data)){
  phenotype.data[i, ncol(phenotype.data)] <- paste(phenotype.data[i,1], phenotype.data[i,2], sep=".")
}

phenotype.data <- phenotype.data[c("SUBJECT.NAME", phenotype.var)]

#drop subjects with any missing phenotype data
phenotype.data <- phenotype.data[!is.na(phenotype.data[,2]),]

#sort phenotype data by subject name
phenotype.data <- phenotype.data[order(phenotype.data[,1]),]

#format dependent variable
y <- phenotype.data[,2]
names(y) <- phenotype.data[,1]

#load diplotype data
diplotype.data <- read.table(cache.file, head=T, colClasses="character")
P <- apply(diplotype.data[,-1], 2, as.numeric)
rownames(P) <- diplotype.data$SUBJECT.NAME

#sort diplotype probabilities by subject name
P <- P[order(rownames(P)),]

#drop subjects without phenotype data
P <- P[rownames(P) %in% names(y),]

#format additive design matrix
A <- matrix(0,nrow=ncol(P),ncol=15)
rownames(A) <- colnames(P)
colnames(A) <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "B1", "B2", "B3", "B4", "B5", "B6", "B7")

for (i in 1:ncol(A)){
  A[i,i] <- 1
}

for (k in (ncol(A)+1):nrow(A)){
  for (i in 1:2){
    A[k,which(colnames(A)==unlist(strsplit(rownames(A)[k],"[.]"))[i])] <- 0.5
  }
}

#specify prior.D
prior.D <- list(P=P, A=A, fixed.diplo=F)

#call TIMBR with full model
#prior.M <- list(model.type="fixed", M.IDs=paste(0:14, collapse=","))
#results.full <- TIMBR(y, prior.D, prior.M, samples=100000)
#save(results.full, file="results_full_locus1.RData")
load("results_full_locus1.RData")

#call TIMBR with CRP model
#prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=1, prior.alpha.rate=2.333415)
#results.crp <- TIMBR(y, prior.D, prior.M, samples=1000000)
#save(results.crp, file="results_crp_locus1.RData")
load("results_crp_locus1.RData")

####################
#report statistics
#calculate median, HPD and length of interval for each approach
library(coda)

colnames(results.full$post.hap.effects) <- colnames(A)
colnames(results.crp$post.hap.effects) <- colnames(A)

HPD.full <- HPDinterval(as.mcmc(results.full$post.hap.effects))
HPD.crp <- HPDinterval(as.mcmc(results.crp$post.hap.effects))

#median
cbind(apply(results.full$post.hap.effects, 2, median), apply(results.crp$post.hap.effects, 2, median))

#HPD
HPD.full
HPD.crp

#HPD length
library(xtable)

cbind(apply(HPD.full, 1, function(x){x[2]-x[1]}), apply(HPD.crp, 1, function(x){x[2]-x[1]}))

table.data <- cbind(apply(HPD.full, 1, function(x){x[2]-x[1]}), apply(HPD.crp, 1, function(x){x[2]-x[1]}))
xtable(table.data, digits=2)

#top allelic series
head(results.crp$p.M.given.y, 10)

table.data <- head(results.crp$p.M.given.y, 10)
table.data <- cbind(sapply(names(table.data), function(x){max(TIMBR:::m.from.M.ID(x))}), table.data)
xtable(table.data, digits=c(0,0,6))

#expectation of number of alleles
mean(results.crp$post.K)

#bayes factor
results.crp$ln.BF - results.full$ln.BF

####################
source("../plot_functions/plot_functions.R")

####################
#generate plots
#combined locus1 plot

pdf(file = "fig_dspr_locus1.pdf", width = 12, height = 3.6)
par(mfrow = c(1, 3))

new_mai <- c(0.6, 0.7, 0.4, 0.15)
par(mai=new_mai)

post.alleles <- rep(0, 15)
names(post.alleles) <- 1:15

for (i in names(results.crp$p.K.given.y)){
  post.alleles[i] <- results.crp$p.K.given.y[i]
}

prior.alleles <- exp(sapply(1:15, TIMBR:::ln.K.prior.crp.marginalized, J=15, a=1, b=2.333415))

df.plot <- barplot(rep(NA,15), 
                   xlab="Number of Alleles", 
                   xaxt='n',
                   cex.lab=1.4, cex.axis=1.2, cex.main=1.5, ylim=c(0,1),
                   las=1)

title(main="Distribution of Alleles", font.main=1, cex.main=1.5)
title(ylab="Probability", line=3.5, cex.lab=1.4)
axis(1, at=df.plot[1 + 0:7*2,1], tick=F, labels=1 + 0:7*2, cex.axis=1.2)

for (h in (1:5)/5){
  abline(h=h, lty=2, col=scales::alpha("black", 0.3))
}

df.plot <- barplot(post.alleles, add=T, 
                   yaxt='n', xaxt='n', ann=FALSE,
                   col=scales::alpha("orange", 0.8))

lines(df.plot, prior.alleles, lty=5, lwd=1.2, col="orange")
points(df.plot, prior.alleles, pch=22, cex=1.2, bg="orange")

legend_homebrew("topright", inset=.02, legend=c("CRP Posterior", "CRP Prior"),
                fill=c(scales::alpha("orange", 0.8), NA),
                border=c("black", "white"),
                bg="white",
                cex=1.2,
                pch=c(NA,22),
                pt.bg=c(NA, "orange"),
                lwd=1.2,
                lty=c(NA,5),
                col=c(NA, "orange"))

box()

put.fig.letter("A", "topleft", font=2, cex=2, offset=c(0.01,-0.01))

hap.names <- colnames(A)
hap.names[8] <- "AB8"

TIMBR.plot.haplotypes(results.full, hap.labels=hap.names, x.lim=c(-5,5), x.lab=phenotype.var,
                      cex.lab=1.4, cex.axis=1.2, cex.main=1.4, font.main = 1,
                      y.lab="",
                      colors=rep(scales::alpha("#4D4D4D", 0.8), 15))
title(main="Full", cex.main=1.5)
title(ylab="Haplotype", line=4, cex.lab=1.4)

box()

put.fig.letter("B", "topleft", font=2, cex=2, offset=c(0.01,-0.01))

TIMBR.plot.haplotypes(results.crp, hap.labels=hap.names, x.lim=c(-5,5), x.lab=phenotype.var,
                      cex.lab=1.4, cex.axis=1.2, cex.main=1.4, font.main = 1,
                      y.lab="",
                      colors=rep(scales::alpha("orange", 0.8),15))
title(main="CRP", cex.main=1.5)
title(ylab="Haplotype", line=4, cex.lab=1.4)

box()

put.fig.letter("C", "topleft", font=2, cex=2, offset=c(0.01,-0.01))

dev.off()
