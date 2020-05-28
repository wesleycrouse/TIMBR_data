phenotype.var <- "CG4086"
cache.file <- "P_locus2.txt"

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
prior.M <- list(model.type="fixed", M.IDs=paste(0:14, collapse=","))
results.full <- TIMBR(y, prior.D, prior.M, samples=100000)
save(results.full, file="results_full_locus2.RData")
#load("results_full_locus2.RData")

#call TIMBR with CRP model
prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=1, prior.alpha.rate=2.333415)
results.crp <- TIMBR(y, prior.D, prior.M, samples=1000000)
save(results.crp, file="results_crp_locus2.RData")
#load("results_crp_locus2.RData")

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
xtable(table.data, digits=c(0,0,4))

#expectation of number of alleles
mean(results.crp$post.K)

#bayes factor
results.crp$ln.BF - results.full$ln.BF

####################
#generate plots
#distribution of number of alleles
post.alleles <- rep(0, 15)
names(post.alleles) <- 1:15

for (i in names(results.crp$p.K.given.y)){
  post.alleles[i] <- results.crp$p.K.given.y[i]
}

prior.alleles <- exp(sapply(1:15, TIMBR:::ln.K.prior.crp.marginalized, J=15, a=1, b=2.333415))

png(filename = "fig_dspr_locus2_alleles.png", width = 600, height = 600)

df.plot <- barplot(post.alleles, xlab="Number of Alleles", ylab="Posterior Probability", main="CG4086", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, ylim=c(0,1))
lines(df.plot, prior.alleles, lty=5)
points(df.plot, prior.alleles, pch=15)

dev.off()

#haplotype effects
hap.names <- colnames(A)
hap.names[8] <- "AB8"

png(filename = "fig_dspr_locus2_haplotypes.png", width = 960, height = 480)

par(mfrow = c(1, 2))

TIMBR.plot.haplotypes(results.full, hap.labels=hap.names, x.lim=c(-4,4), x.lab=phenotype.var)
title(main="Full Model")

TIMBR.plot.haplotypes(results.crp, hap.labels=hap.names, x.lim=c(-4,4), x.lab=phenotype.var)
title(main="CRP Model")

dev.off()
