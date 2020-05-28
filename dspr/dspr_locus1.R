##### identify peak loci

lod <- read.table("data/LODscores_eqtl_trimmed.txt", head=T)
lod[which.max(lod[,3]),]
lod[which.max(lod[,4]),]

##### load diplotype probabilities for each RIL

P.A <- read.table("data/HMMregA_R2_trimmed1.txt")
P.B <- read.table("data/HMMregB_R2_trimmed1.txt")

P.A <- P.A[,1:39]
P.B <- P.B[,1:39]

colnames(P.A) <- c("Chr", "Pos", "RIL.A", "A1.A1", "A1.A2", "A1.A3","A1.A4","A1.A5","A1.A6","A1.A7","A1.A8",
  "A2.A2","A2.A3","A2.A4","A2.A5","A2.A6","A2.A7","A2.A8",
  "A3.A3","A3.A4","A3.A5","A3.A6","A3.A7","A3.A8",
  "A4.A4","A4.A5","A4.A6","A4.A7","A4.A8",
  "A5.A5","A5.A6","A5.A7","A5.A8",
  "A6.A6","A6.A7","A6.A8",
  "A7.A7","A7.A8",
  "A8.A8")

colnames(P.B) <- c("Chr", "Pos", "RIL.B", "B1.B1", "B1.B2", "B1.B3","B1.B4","B1.B5","B1.B6","B1.B7","B1.A8",
  "B2.B2","B2.B3","B2.B4","B2.B5","B2.B6","B2.B7","B2.A8",
  "B3.B3","B3.B4","B3.B5","B3.B6","B3.B7","B3.A8",
  "B4.B4","B4.B5","B4.B6","B4.B7","B4.A8",
  "B5.B5","B5.B6","B5.B7","B5.A8",
  "B6.B6","B6.B7","B6.A8",
  "B7.B7","B7.A8",
  "A8.A8")

#####load phenotype file to identify crosses for RIX
pheno <- read.table("data/FemaleHeadExpression_trimmed.txt", head=T)

#####construct diplotype probability matrix for RIX
P <- matrix(0, nrow(pheno), 15 + choose(15,2))

founders <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "B1", "B2", "B3", "B4", "B5", "B6", "B7")
colnames(P) <- rep(NA, 15 + choose(15,2))

for (i in 1:15){
  colnames(P)[i] <- paste(founders[i], founders[i], sep=".")
}

k <- 15
for (i in 1:14){
  for (j in (i+1):15){
    k <- k+1
    colnames(P)[k] <- paste(founders[i], founders[j], sep=".")
  }
}

#calculate diplotype probabilities
for (n in 1:nrow(P)){
  print(n)
  
  a <- which(pheno[n,2]==P.A[,3])
  b <- which(pheno[n,1]==P.B[,3])
  
  for (i in 4:ncol(P.A)){
    for (j in 4:ncol(P.B)){
      for (k in 1:2){
        for (l in 1:2){
          current.diplo <- paste(unlist(strsplit(colnames(P.A)[i], "[.]"))[k], unlist(strsplit(colnames(P.B)[j], "[.]"))[l], sep=".")
          P[n,colnames(P)==current.diplo] <- P[n,colnames(P)==current.diplo] + 0.5*P.A[a,i]*0.5*P.B[b,j]
        }
      }
    }
  }
}

for (n in 1:nrow(P)){
  P[n,] <- P[n,]/sum(P[n,])
}

SUBJECT.NAME <- rep(NA, nrow(P))

for (i in 1:nrow(P)){
  SUBJECT.NAME[i] <- paste(pheno[i,1], pheno[i,2], sep=".")
}

P <- cbind(SUBJECT.NAME, P)

write.table(P, "P_locus1.txt", col.names=T, quote=F)
