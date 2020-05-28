load("TIMBR_sim_results_full.RData")

library(ggplot2)
library(TIMBR)
library(gridExtra)

####################
#function for binomial confidence interval
binomial.prop.ci <- function(p, n.sims=1000, alpha=0.05){
  ci.lower <- qbeta(0.5*alpha, p*n.sims+0.5, (1-p)*n.sims+0.5)
  ci <- c(ci.lower, qbeta(1-0.5*alpha, p*n.sims+0.5, (1-p)*n.sims+0.5))
  
  if (p==0){
    ci[1] <- 0
  } else if (p==1){
    ci[2] <- 1
  }
  
  ci
}

####################
#prior number of alleles

hyperparam <- calc.concentration.prior(8, 0.05, 0.01)
prior.gamma <- exp(sapply(1:8, TIMBR:::ln.K.prior.crp.marginalized, J=8, a=hyperparam[1], b=hyperparam[2]))
prior.exp <- exp(sapply(1:8, TIMBR:::ln.K.prior.crp.marginalized, J=8, a=1, b=2.333415))

prior.uniform <- table(sapply(ewenss.calc(8, list(type="gamma", shape=1, rate=1))$M.IDs, 
                              function(x){max(as.numeric(unlist(strsplit(x, split=","))))+1}))/exp(TIMBR:::ln.bell(8))

png(filename = "fig_3_0_prior.png", width = 480, height = 480)

plot(1:8, prior.gamma, type="l", lty=5, col="grey", ylim=c(0,1), ylab="Prior Probability", xlab="Number of Alleles", main="Prior Distribution of Alleles")
points(1:8, prior.gamma, pch=22, bg="grey")
lines(1:8, prior.exp, lty=5, col="orange")
points(1:8, prior.exp, pch=22, bg="orange")
lines(1:8, prior.uniform, lty=5, col="darkcyan")
points(1:8, prior.uniform, pch=22, bg="darkcyan")

legend(1, legend=c("Uniform", "Gamma", "Exponential"), col=c("darkcyan", "grey", "orange"), lty=5)

dev.off()

####################
#accuracy - (0-1)

colors <- c("grey", "orange", "darkcyan")
strategies <- c("gamma.st.05.01", "exponential.st.5", "uniform")
accuracy <- aggregate(accuracy ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=mean)

png(filename = "fig_3_1_accuracy.png", width = 960, height = 480)

par(mfrow = c(1, 2))

for (v in c(0.1, 0.5)){
  for (A in c(1)){
    plot(c(), c(), ylim=c(0,1), xlim=c(1,8), las=1, xlab="Number of Alleles", ylab = "0-1 Accuracy", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
         main=paste0("QTL Effect Size: ", v))
    
    for (h in (1:4)/5){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3))
    }
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "accuracy")]
      data.subset <- cbind(data.subset, t(sapply(data.subset[,2], binomial.prop.ci)))
      colnames(data.subset)[3:4] <- c("upper", "lower")
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
      
    }
  }
}

legend(1, 0.2, legend=c("Uniform", "Gamma", "Exponential"), col=c("darkcyan", "grey", "orange"), lty=1)

dev.off()

####################
#accuracy - posterior.M

accuracy <- aggregate(posterior.M ~ strategy + alleles + alpha + var.exp, data=results.full[results.full$strategy %in% strategies,], FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
accuracy <- cbind(accuracy[,1:4], accuracy$posterior.M)
colnames(accuracy)[-c(1:4)] <- c("posterior.M", "lower", "upper")

png(filename = "fig_3_2_posterior_M.png", width = 960, height = 480)

par(mfrow = c(1, 2))

for (v in c(0.1, 0.5)){
  for (A in c(1)){
    plot(c(), c(), ylim=c(0,1), xlim=c(1,8), las=1, xlab="Number of Alleles", ylab = "Posterior Certainty", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
         main=paste0("QTL Effect Size: ", v))
    
    for (h in (1:4)/5){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3))
    }
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "posterior.M", "lower", "upper")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
      
    }
  }
}

legend(6, 1, legend=c("Uniform", "Gamma", "Exponential"), col=c("darkcyan", "grey", "orange"), lty=1)

dev.off()

####################
#accuracy - exp.K

accuracy <- aggregate(exp.K ~ strategy + alleles + alpha + var.exp, data=results.full[results.full$strategy %in% strategies,], FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
accuracy <- cbind(accuracy[,1:4], accuracy$exp.K)
colnames(accuracy)[-c(1:4)] <- c("exp.K", "lower", "upper")

png(filename = "fig_3_3_exp_K.png", width = 960, height = 480)

par(mfrow = c(1, 2))

for (v in c(0.1, 0.5)){
  for (A in c(1)){
    plot(c(), c(), ylim=c(1,8), xlim=c(1,8), las=1, xlab="Number of Alleles", ylab = "Expected Number of Alleles", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
         main=paste0("QTL Effect Size: ", v))
    
    abline(0, 1, lty=3)
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "exp.K", "upper", "lower")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
    }
  }
}

legend(5.5, 2.5, legend=c("Uniform", "Gamma", "Exponential"), col=c("darkcyan", "grey", "orange"), lty=1)

dev.off()

####################
#error - MSE

colors <- c("grey", "orange", "darkcyan", "green", "black")
strategies <- c("gamma.st.05.01", "exponential.st.5", "uniform", "known.M", "full")

error <- aggregate(MSE ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
error <- cbind(error[,1:4], error$MSE)
colnames(error)[-c(1:4)] <- c("MSE", "lower", "upper")

png(filename = "fig_3_5_MSE.png", width = 960, height = 480)

par(mfrow = c(1, 2))

for (v in c(0.1, 0.5)){
  for (A in c(1)){
    plot(c(), c(), ylim=c(0,0.4), xlim=c(1, 8), las=1, xlab="Number of Alleles", ylab = "Mean Squared Error", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
         main=paste0("QTL Effect Size: ", v))
    
    for (h in (1:4)/10){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3))
    }
    
    for (i in 1:length(strategies)){
      data.subset <- error[error$strategy==strategies[i] & error$alpha==A & error$var.exp==v, c("alleles", "MSE", "upper", "lower")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
    }
  }
}

legend(1, 0.4, legend=c("Uniform", "Gamma", "Exponential", "Full", "Known"), col=c("darkcyan", "grey", "orange", "black", "green"), lty=1)

dev.off()

####################
#accuracy - (0-1) - exponential

colors <- c("orange", "blue", "red", "purple")
strategies <- c("exponential.st.5", "exponential.st.5.tree", "exponential.st.5.tree.wrong", "exponential.st.5.tree.miss")
accuracy <- aggregate(accuracy ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=mean)

A <- 1

png(filename = paste0("fig_3_1_e_", A, "_accuracy.png"), width = 960, height = 480,)

par(mfrow = c(1, 2))

v <- 0.1

plot(c(), c(), ylim=c(0,1), xlim=c(1,8), las=1, xlab="Number of Alleles", ylab = "0-1 Accuracy", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
     main=paste0("QTL Effect Size: ", v))

for (h in (1:4)/5){
  abline(h=h, lty=2, col=scales::alpha("black", 0.3))
}

for (i in 1:length(strategies)){
  data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "accuracy")]
  data.subset <- cbind(data.subset, t(sapply(data.subset[,2], binomial.prop.ci)))
  colnames(data.subset)[3:4] <- c("upper", "lower")
  
  polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
  lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
  points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
}

legend(5.5, 1, legend=c("CRP", "Tree", "Misspecified", "Incorrect"), col=c("orange", "blue", "purple", "red"), lty=1)

v <- 0.5

plot(c(), c(), ylim=c(0,1), xlim=c(1,8), las=1, xlab="Number of Alleles", ylab = "0-1 Accuracy", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
     main=paste0("QTL Effect Size: ", v))

for (h in (1:4)/5){
  abline(h=h, lty=2, col=scales::alpha("black", 0.3))
}

for (i in 1:length(strategies)){
  data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "accuracy")]
  data.subset <- cbind(data.subset, t(sapply(data.subset[,2], binomial.prop.ci)))
  colnames(data.subset)[3:4] <- c("upper", "lower")
  
  polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
  lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
  points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
}

dev.off()

####################
#error - MSE - exponential

colors <- c("orange", "blue", "red", "purple", "green", "black")
strategies <- c("exponential.st.5", "exponential.st.5.tree", "exponential.st.5.tree.wrong", "exponential.st.5.tree.miss", "known.M", "full")

error <- aggregate(MSE ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
error <- cbind(error[,1:4], error$MSE)
colnames(error)[-c(1:4)] <- c("MSE", "lower", "upper")

A <- 1

png(filename = paste0("fig_3_5_e_", A, "_MSE.png"), width = 960, height = 480)

par(mfrow = c(1, 2))

for (v in c(0.1, 0.5)){
  plot(c(), c(), ylim=c(0,0.4), xlim=c(1, 8), las=1, xlab="Number of Alleles", ylab = "Mean Squared Error", frame.plot=T, cex.lab=1.5, cex.axis=1.3,
       main=paste0("QTL Effect Size: ", v))
  
  for (h in (1:4)/10){
    abline(h=h, lty=2, col=scales::alpha("black", 0.3))
  }
  
  for (i in 1:length(strategies)){
    data.subset <- error[error$strategy==strategies[i] & error$alpha==A & error$var.exp==v, c("alleles", "MSE", "upper", "lower")]
    
    polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
    lines(data.subset$alleles, data.subset[,2], lwd=3, col=colors[i])
    points(data.subset$alleles, data.subset[,2], col=colors[i], pch=16)
  }
}

legend(1, 0.4, legend=c("CRP", "Tree", "Misspecified", "Incorrect", "Full", "Known"), col=c("orange", "blue", "purple", "red", "black", "green"), lty=1)

dev.off()
