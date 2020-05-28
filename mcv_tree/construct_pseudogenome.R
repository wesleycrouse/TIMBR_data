#specify command line options and parse command line arguments

options(echo=FALSE)
options(stringsAsFactors = FALSE)

command.line <- commandArgs(trailingOnly = TRUE)

job.id <- command.line[1]
chr <- command.line[2]
start <- command.line[3]
end <- command.line[4]

#specify start and end positions
#job.id <- "chr7_3776"
#chr <- "7"
#start <- "103807679"
#end <- "103831178"

#transform start and end positions from strings
start <- as.numeric(start)
end <- as.numeric(end)

#load libraries

#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("VariantAnnotation")

library(Biostrings)
library(VariantAnnotation)
library(stringr)

#specify filename of VCF files and reference
snp.vcf.file <- paste(job.id, "_snps.recode.vcf", sep="")
indel.vcf.file <- paste(job.id, "_indels.recode.vcf", sep="")
reference.file <- paste("Mus_musculus.GRCm38.dna.chromosome.", chr, ".fa", sep="")

#load reference sequence
ref <- readDNAStringSet(reference.file, format="fasta")
ref <- DNAStringSet(ref, start=start, end=end)
names(ref) <- "ref"

#load SNP genotypes as nucleotides
GT.snps <- readGT(snp.vcf.file, nucleotides=TRUE)
GT.snps <- matrix(sapply(strsplit(GT.snps,"/"), "[", 1), nrow(GT.snps), ncol(GT.snps))

#load SNP VCF file
snps <- readVcf(snp.vcf.file, genome="mm10")

#store SNP quality indicator and set quality indicator to 0 if equal to reference
FI.snps <- geno(snps)$FI
FI.snps[is.na(FI.snps)] <- 0

FI.snps <- FI.snps
FI.snps[geno(snps)$GT=="0/0"] <- 0

#store SNP sample names
strains.snps <- samples(header(snps))

#store SNP locations and adjust for start position
pos.snps <- IRanges(snps@rowRanges@ranges@start-start+1, width=snps@rowRanges@ranges@width)
#pos.snps <- IRanges(snps@rowData@ranges@start-start+1, width=snps@rowData@ranges@width)

#copy reference for each pseudogenome
pseudogenomes <- rep(ref, length(strains.snps))
names(pseudogenomes) <- strains.snps

#replace SNPs that satisfy the quality filter for each pseudogenome
for (i in 1:length(strains.snps)){
  pseudogenomes[i] <- replaceAt(pseudogenomes[i], at=pos.snps[FI.snps[,i]==1], GT.snps[FI.snps[,i]==1,i])
}

###################

#load indel VCF file
indels <- readVcf(indel.vcf.file, genome="mm10")

if (nrow(geno(indels)$GT)!=0){
  #load indel genotypes as nucleotides
  GT.indels <- readGT(indel.vcf.file, nucleotides=TRUE)
  GT.indels <- matrix(sapply(strsplit(GT.indels,"/"), "[", 1), nrow(GT.indels), ncol(GT.indels))
  
  #store indel quality indicator
  FI.indels <- geno(indels)$FI
  FI.indels[is.na(FI.indels)] <- 0
  
  #store indel sample names and ensure that pseudogenomes are in the correct order
  strains.indels <- samples(header(indels))
  pseudogenomes <- pseudogenomes[order(strains.indels)]
  
  #store indel locations and adjust for start position
  pos.indels <- IRanges(indels@rowRanges@ranges@start-start+1, width=indels@rowRanges@ranges@width)
  #pos.indels <- IRanges(indels@rowData@ranges@start-start+1, width=indels@rowData@ranges@width)

  #filter indels and adjust lengths
  flag <- rep(0, nrow(GT.indels))
  
  for (i in 1:nrow(GT.indels)){
    #snp in region of indel
    if (length(intersect(pos.indels[i],pos.snps[rowSums(FI.snps)!=0]))>0){
      flag[i] <- 4
    }
    #complex deletion
    if (ref(indels)@ranges@width[i]>1 & max(c(unlist(sapply(GT.indels[i, FI.indels[i,]!=0 & GT.indels[i,]!=ref(indels)[i]], nchar)),1))!=1){
      flag[i] <- 3
    }
    #multiple alternate alleles
    if (nlevels(as.factor(GT.indels[i,GT.indels[i,]!=ref(indels)[i] & FI.indels[i,]==1]))>1){
      flag[i] <- 2
    }
    #no high-confidence non-reference alleles
    if (sum(FI.indels[i,GT.indels[i,]!=ref(indels)[i]])==0){
      flag[i] <- 1
    }
    #process insertions and deletions
    if (ref(indels)@ranges@width[i]>1 & flag[i]==0){
      GT.indels[i,FI.indels[i,]==0] <- as.character(ref(indels)[i])
      GT.indels[i,FI.indels[i,]!=0] <- str_pad(GT.indels[i,FI.indels[i,]!=0], ref(indels)@ranges@width[i], side="right", pad="-")
    } else if (flag[i]==0) {
      same.as.ref <- GT.indels[i,]==ref(indels)[i]
      current.length <- max(sapply(GT.indels[i, FI.indels[i,]!=0], nchar))
      ref(indels)[i] <- str_pad(as.character(ref(indels)[i]), current.length, side="right", pad="-")
      GT.indels[i, FI.indels[i,]==0 | same.as.ref] <- as.character(ref(indels)[i])
    }
  }
  
  #discard indels that overlap keeping the first one
  flag[flag==0][disjointBins(pos.indels[flag==0])!=1] <- 5

  #discard indels that span endpoint
  flag[pos.indels@start+pos.indels@width-1 > ref@ranges@width] <- 6

  #replace indels that satisfy the filter criteria for each pseudogenome and adjust reference sequence
  for (i in 1:length(strains.indels)){
    pseudogenomes[i] <- replaceAt(pseudogenomes[i], at=pos.indels[flag==0], GT.indels[flag==0,i])
  }
  
  ref <- replaceAt(ref, at=pos.indels[flag==0], ref(indels)[flag==0])
}

#rename sequences based on CC founder abbreviations, add reference to object and reorder
names(pseudogenomes)[names(pseudogenomes)=="A_J"] <- "A"
names(ref) <- "B"
names(pseudogenomes)[names(pseudogenomes)=="129S1_SvImJ"] <- "C"
names(pseudogenomes)[names(pseudogenomes)=="NOD_ShiLtJ"] <- "D"
names(pseudogenomes)[names(pseudogenomes)=="NZO_HlLtJ"] <- "E"
names(pseudogenomes)[names(pseudogenomes)=="CAST_EiJ"] <- "F"
names(pseudogenomes)[names(pseudogenomes)=="PWK_PhJ"] <- "G"
names(pseudogenomes)[names(pseudogenomes)=="WSB_EiJ"] <- "H"

pseudogenomes <- c(pseudogenomes, ref)
pseudogenomes <- pseudogenomes[order(names(pseudogenomes))]

#write FASTA
writeXStringSet(pseudogenomes, paste(job.id, ".fa", sep=""), append=FALSE)

message(paste(sum(rowSums(FI.snps)!=0), " SNPs used", sep=""))
if (nrow(geno(indels)$GT)!=0){
  message(paste(sum(flag==0), " indels used", sep=""))
} else {
  message("0 indels used")
}
