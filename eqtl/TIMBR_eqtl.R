command.line <- commandArgs(trailingOnly = TRUE)

phenotype.var <- as.character(command.line[1]) 
locus <- as.character(command.line[2])

#phenotype.var <- "X1212612"
#locus <- "M17.345"

#drop \r that appears when running from command line
locus <- unlist(strsplit(locus, "[^[:print:]]"))

####################
library(TIMBR)
library(data.table)

#load phenotype data
phenotype.data <- fread("data/rankz-filtered-precc-geneexp.txt", header=T, stringsAsFactors=F, select=c("SUBJECT.NAME", phenotype.var))
phenotype.data <- as.data.frame(phenotype.data)

#drop subjects with any missing phenotype data
phenotype.data <- phenotype.data[!is.na(phenotype.data[,2]),]

#sort phenotype data by subject name
phenotype.data <- phenotype.data[order(phenotype.data[,1]),]

#format dependent variable
y <- phenotype.data[,2]
names(y) <- phenotype.data[,1]

#load diplotype data
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

#call TIMBR with CRP model
prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=1, prior.alpha.rate=2.333415)
results.crp <- TIMBR(y, prior.D, prior.M, samples=100000)

save(results.crp, file=paste0("results/results_crp_", phenotype.var, ".RData"))
