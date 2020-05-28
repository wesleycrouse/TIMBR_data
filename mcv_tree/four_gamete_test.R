#specify command line options and parse command line arguments
options(echo=FALSE)
options(stringsAsFactors = FALSE)

command.line <- commandArgs(trailingOnly = TRUE)

chr <- command.line[1]

#set job id
job.id <- paste("chr", chr, sep="")

#load libraries

#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("VariantAnnotation")

library(Biostrings)
library(VariantAnnotation)
library(stringr)

#specify filename of VCF files and reference
snp.vcf.file <- paste(job.id, "_snps.recode.vcf", sep="")

#load portions of VCF file
param <- ScanVcfParam(geno=c("GT","FI"))
snps <- readVcf(snp.vcf.file, genome="mm10", param=param)

#load SNP genotypes as nucleotides
GT.snps <- geno(snps)$GT
GT.snps <- matrix(sapply(strsplit(GT.snps,"/"), "[", 1), nrow(GT.snps), ncol(GT.snps))

#store SNP quality indicator and set quality indicator to 0 if equal to reference
FI.snps <- geno(snps)$FI
FI.snps[is.na(FI.snps)] <- 0

FI.snps <- FI.snps
FI.snps[geno(snps)$GT=="0/0"] <- 0

#set genotype equal to reference if low quality
for (i in 1:nrow(GT.snps)){
  GT.snps[i,FI.snps[i,]==0] <- 0
}

#store SNP locations and adjust for start position
pos.snps <- IRanges(snps@rowData@ranges@start, width=snps@rowData@ranges@width)

#remove unused information from memory
rm(snps)

#breakpoints
breakpoints <- c()

#perform four-gamete test
for (i in 1:(nrow(GT.snps)-1)){
  
  if (i %% 1000 == 0){
    print(i)
  }
  
  #ensure that both snps are biallelic
  if (nlevels(as.factor(c(GT.snps[i,],"0")))<=2 & nlevels(as.factor(c(GT.snps[i+1,],"0")))<=2){
    haplotypes <- as.factor(apply(rbind(c(GT.snps[i,],"0"), c(GT.snps[i+1,],"0")),2,paste, collapse=""))
    if (nlevels(haplotypes)==4){
      breakpoints <- c(breakpoints, i)
    }
  }
}

#set start and end by first and last SNPs
start <- pos.snps[1,]@start
end <- pos.snps[length(pos.snps),]@start

#calculate midpoints between SNPs that demonstrate evidence of recombination
end.breakpoints <- pos.snps[breakpoints]@start + floor((pos.snps[breakpoints+1]@start - pos.snps[breakpoints]@start)/2)
end.breakpoints <- c(end.breakpoints, end)

start.breakpoints <- pos.snps[breakpoints]@start + floor((pos.snps[breakpoints+1]@start - pos.snps[breakpoints]@start)/2)
start.breakpoints <- c(start, start.breakpoints+1)

tree.ranges <- cbind(start.breakpoints, end.breakpoints)

#store breakpoints
#write.table(tree.ranges, paste(job.id, "_ranges.txt", sep=""), quote=FALSE, row.names=FALSE)
