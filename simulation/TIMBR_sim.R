args <- commandArgs(TRUE)

job.id <- as.character(args[1])
N.J <- as.numeric(args[2])
K <- as.numeric(args[3])
true.alpha <- as.numeric(args[4])
var.exp <- as.numeric(args[5])
iter <- as.numeric(args[6])

#job.id <- "test"
#N.J <- 50
#K <- 3
#true.alpha <- 1
#var.exp <- 0.5
#iter <- 2

####################
#libraries and constant settings for all simulations
#devtools::install_github("wesleycrouse/TIMBR")
library(TIMBR)
library(ape)

J <- 8
samples <- 100000

##################
#create object and function to store results
results.df <- data.frame(iter=as.numeric(), strategy=as.character(), accuracy=logical(), MSE=numeric(), ln.BF=numeric(), exp.K=numeric(), posterior.M=numeric(), prob.over=numeric(), prob.under=numeric(), stringsAsFactors=F)

update.results <- function(strategy, results){
  this.result <- list(iter=i, strategy=strategy, accuracy=NA, MSE=NA, ln.BF=NA, exp.K=NA, posterior.M=NA, prob.over=NA, prob.under=NA)
  
  this.result$accuracy <- ifelse(names(results$p.M.given.y[1])==true.M.ID, 1, 0)
  this.result$MSE <- mean(colSums((t(results$post.hap.effects) - true.B)^2))
  this.result$ln.BF <- results$ln.BF
  this.result$exp.K <- sum(as.numeric(names(results$p.K.given.y))*results$p.K.given.y)
  this.result$posterior.M <- results$p.M.given.y[true.M.ID]
  this.result$prob.over <- mean(results$post.K>K)
  this.result$prob.under <- mean(results$post.K<K)
  
  results.df[nrow(results.df)+1,] <<- this.result
}

##################
#precomputed values
hyperparam <- calc.concentration.prior(8, 0.05, 0.01)

##################

print(job.id)

for (i in 1:iter){
  print(i)
  
  #simulate data
  tree <- rcoal(J, LETTERS[1:J])
  
  #save(tree, file=paste("results/", job.id, "_tree.RData", sep=""))
  
  true.prior.M <- ewenss.calc(tree, list(type="fixed", alpha=true.alpha))
  
  which.K <- sapply(true.prior.M$M.IDs, function(x){max(TIMBR:::m.from.M.ID(x))})==K
  M.prob <- true.prior.M$ln.probs[which.K]
  names(M.prob) <- true.prior.M$M.IDs[which.K]
  M.prob <- exp(M.prob - matrixStats::logSumExp(M.prob))
  
  true.M.ID <- sample(names(M.prob), 1, prob=M.prob)
  true.M <- TIMBR:::M.matrix.from.ID(true.M.ID)
  
  sim.data <- TIMBR:::simulate.population(true.M, N.J, var.exp)
  
  true.B <- c(true.M%*%sim.data$B)
  
  y <- sim.data$y
  prior.D <- sim.data$prior.D
  
  rm(M.prob, which.K)
  
  ##################
  #simulate incorrect trees
  #same tree with total distance resampled from coalesecent
  #tree.dist <- tree
  #tree.dist$edge.length <- tree.dist$edge.length/sum(tree.dist$edge.length)*sum(rcoal(J, LETTERS[1:J])$edge.length)
  
  #scrambled tree based on Ansari 2016 with total distance resampled independently from coalesecent
  tree.miss <- cophenetic(tree)
  tree.miss[upper.tri(tree.miss)] <- sapply(tree.miss[upper.tri(tree.miss)], function(x){runif(1, x*0.1, x*0.9)})
  tree.miss[lower.tri(tree.miss)] <- t(tree.miss)[lower.tri(tree.miss)]
  tree.miss <- as.phylo(hclust(as.dist(tree.miss), "average"))
  
  tree.miss$edge.length <- tree.miss$edge.length/sum(tree.miss$edge.length)*sum(rcoal(J, LETTERS[1:J])$edge.length)
  
  #unrelated tree sampled from coalesecent
  tree.wrong <- rcoal(J, LETTERS[1:J])
  
  ##################
  #analysis using TIMBR
  
  ##################
  strategy <- "known.M"
  print(strategy)
  prior.M <- list(model.type="fixed", M.IDs=true.M.ID)
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "full"
  print(strategy)
  prior.M <- list(model.type="fixed", M.IDs=paste(0:(J-1), collapse=","))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "gamma.st.05.01"
  print(strategy)
  prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=hyperparam[1], prior.alpha.rate=hyperparam[2])
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "uniform"
  print(strategy)
  prior.M <- list(model.type="uniform")
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "exponential.st.5"
  print(strategy)
  prior.M <- list(model.type="crp", prior.alpha.type="gamma", prior.alpha.shape=1, prior.alpha.rate=2.333415)
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "gamma.st.05.01.tree"
  print(strategy)
  prior.M <- ewenss.calc(tree, list(type="gamma", shape=hyperparam[1], rate=hyperparam[2]))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)

  ##################
  strategy <- "gamma.st.05.01.tree.wrong"
  print(strategy)
  prior.M <- ewenss.calc(tree.wrong, list(type="gamma", shape=hyperparam[1], rate=hyperparam[2]))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)

  ##################
  strategy <- "gamma.st.05.01.tree.miss"
  print(strategy)
  prior.M <- ewenss.calc(tree.miss, list(type="gamma", shape=hyperparam[1], rate=hyperparam[2]))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "exponential.st.5.tree"
  print(strategy)
  prior.M <- ewenss.calc(tree, list(type="gamma", shape=1, rate=2.333415))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "exponential.st.5.tree.wrong"
  print(strategy)
  prior.M <- ewenss.calc(tree.wrong, list(type="gamma", shape=1, rate=2.333415))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
  
  ##################
  strategy <- "exponential.st.5.tree.miss"
  print(strategy)
  prior.M <- ewenss.calc(tree.miss, list(type="gamma", shape=1, rate=2.333415))
  results <- TIMBR(y, prior.D, prior.M, verbose=F, samples = samples)
  update.results(strategy, results)
}

save(results.df, file=paste("results/", job.id, ".RData", sep=""))
