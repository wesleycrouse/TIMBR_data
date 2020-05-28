NJ <- 50
alleles <- 1:8
jobs <- 1:100
alpha <- c(0.1, 1, 10)
var.exp <- c(0.1, 0.5)

for (K in alleles){
  for (A in alpha){
    for (v in var.exp){
      for (i in jobs){
        job.id <- paste0("results/TIMBR_sim_N.J_", NJ, "_K_", K, "_alpha_", A, "_v_", v, "_job_", i, ".RData")
        load(job.id)
        
        #coerce to matrix because of unexplained error when rbinding as data frames
        results.df <- as.matrix(cbind(results.df, K, A, v, i))

        if (!exists("results.full")){
          results.full <- results.df
        } else {
          results.full <- rbind(results.full, results.df)
        }
      }
    }
  }
}

colnames(results.full)[(length(colnames(results.full))-3):length(colnames(results.full))] <- c("alleles", "alpha", "var.exp", "job")

#coerce back to data frame and set variable classes
results.full <- as.data.frame(results.full, stringsAsFactors=F)

results.full$iter <- as.integer(results.full$iter)
results.full$accuracy <- as.numeric(results.full$accuracy)
results.full$MSE <- as.numeric(results.full$MSE)
results.full$ln.BF <- as.numeric(results.full$ln.BF)
results.full$exp.K <- as.numeric(results.full$exp.K)
results.full$posterior.M <- as.numeric(results.full$posterior.M)
results.full$prob.over <- as.numeric(results.full$prob.over)
results.full$prob.under <- as.numeric(results.full$prob.under)
results.full$alleles <- as.integer(results.full$alleles)
results.full$alpha<- as.numeric(results.full$alpha)
results.full$var.exp <- as.numeric(results.full$var.exp)
results.full$job <- as.integer(results.full$job)

save(results.full, file="TIMBR_sim_results_full.RData")
