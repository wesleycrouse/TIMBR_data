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

#function for legend
source("../plot_functions/plot_functions.R")

####################
#prior number of alleles

hyperparam <- calc.concentration.prior(8, 0.05, 0.01)
prior.gamma <- exp(sapply(1:8, TIMBR:::ln.K.prior.crp.marginalized, J=8, a=hyperparam[1], b=hyperparam[2]))
prior.exp <- exp(sapply(1:8, TIMBR:::ln.K.prior.crp.marginalized, J=8, a=1, b=2.333415))

prior.uniform <- table(sapply(ewenss.calc(8, list(type="gamma", shape=1, rate=1))$M.IDs, 
                              function(x){max(as.numeric(unlist(strsplit(x, split=","))))+1}))/exp(TIMBR:::ln.bell(8))

pdf(file = "fig_3_0_prior.pdf", width=4, height=4)

par(mai=c(0.75, 0.8, 0.5, 0.2))

plot(c(), c(), 
     ylim=c(0,1), 
     xlim=c(1,8),
     ylab="Prior Probability", 
     xlab="",
     main="Distribution of Alleles", 
     font.main=1,
     cex.lab=1.1, 
     cex.axis=1, 
     cex.main=1.2, 
     lwd=1, 
     las=1)

title(xlab="Number of Alleles", line=2.5, cex.lab=1.1)

for (h in (1:5)/5){
  abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
}

lines(1:8, prior.gamma, lty=5, lwd=1, col="grey")
points(1:8, prior.gamma, pch=22, cex=1, bg="grey")
lines(1:8, prior.exp, lty=5, lwd=1, col="orange")
points(1:8, prior.exp, pch=22, cex=1, bg="orange")
lines(1:8, prior.uniform, lty=5, lwd=1, col="darkcyan")
points(1:8, prior.uniform, pch=22, cex=1, bg="darkcyan")

legend_homebrew("topright", inset=.02, legend=c("Uniform", "Gamma", "Exponential"),
                bg="white",
                cex=0.8,
                pch=22,
                pt.bg=c("darkcyan", "grey", "orange"),
                lwd=1,
                lty=5,
                col=c("darkcyan", "grey", "orange"))

dev.off()

####################
#accuracy - (0-1)

colors <- c("grey", "orange", "darkcyan")
strategies <- c("gamma.st.05.01", "exponential.st.5", "uniform")
accuracy <- aggregate(accuracy ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=mean)

pdf(file = "fig_3_1_accuracy.pdf", width = 8, height = 4)

par(mfrow = c(1, 2))
par(mai=c(0.8, 0.8, 0.5, 0.2))

counter <- 0

for (v in c(0.1, 0.5)){
  for (A in 1){
    counter <- counter + 1
    
    plot(c(), c(), 
         ylim=c(0,1), 
         xlim=c(1,8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1.1, 
         font.main = 1, 
         frame.plot=T,
         xlab="", ylab = "0-1 Accuracy",
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    for (h in (1:5)/5){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
    }
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "accuracy")]
      data.subset <- cbind(data.subset, t(sapply(data.subset[,2], binomial.prop.ci)))
      colnames(data.subset)[3:4] <- c("upper", "lower")
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
  }
  
  if (v==0.1){
    legend("topright", inset=.02, legend=c("Uniform", "Gamma", "Exponential"),
           bg="white",
           cex=0.8,
           pch=21,
           pt.bg=c("darkcyan", "grey", "orange"),
           lwd=1.5,
           lty=1,
           col=c("darkcyan", "grey", "orange"))
  }
}

dev.off()

####################
#accuracy - posterior.M

accuracy <- aggregate(posterior.M ~ strategy + alleles + alpha + var.exp, data=results.full[results.full$strategy %in% strategies,], FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
accuracy <- cbind(accuracy[,1:4], accuracy$posterior.M)
colnames(accuracy)[-c(1:4)] <- c("posterior.M", "lower", "upper")

pdf(file = "fig_3_2_posterior_M.pdf", width = 8, height = 4)

par(mfrow = c(1, 2))
par(mai=c(0.8, 0.8, 0.5, 0.2))

counter <- 0

for (v in c(0.1, 0.5)){
  for (A in 1){
    counter <- counter + 1
    
    plot(c(), c(), 
         ylim=c(0,1), 
         xlim=c(1,8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1.1, 
         font.main = 1, 
         frame.plot=T,
         xlab="", 
         ylab = "Posterior Certainty", 
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    for (h in (1:5)/5){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
    }
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "posterior.M", "lower", "upper")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
  }
}

legend("topright", inset=.02, legend=c("Uniform", "Gamma", "Exponential"),
                bg="white",
                cex=0.8,
                pch=21,
                pt.bg=c("darkcyan", "grey", "orange"),
                lwd=1.5,
                lty=1,
                col=c("darkcyan", "grey", "orange"))

dev.off()

####################
#accuracy - exp.K

accuracy <- aggregate(exp.K ~ strategy + alleles + alpha + var.exp, data=results.full[results.full$strategy %in% strategies,], FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
accuracy <- cbind(accuracy[,1:4], accuracy$exp.K)
colnames(accuracy)[-c(1:4)] <- c("exp.K", "lower", "upper")

pdf(file = "fig_3_3_exp_K.pdf", width = 8, height = 4)

par(mfrow = c(1, 2))
par(mai=c(0.8, 0.8, 0.5, 0.2))

counter <- 0

for (v in c(0.1, 0.5)){
  for (A in 1){
    counter <- counter + 1
    
    plot(c(), c(), ylim=c(1,8), xlim=c(1,8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1.1, 
         font.main = 1, 
         frame.plot=T,
         xlab="", 
         ylab = "Expected Number of Alleles",
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    abline(0, 1, lty=3, lwd=0.7)
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "exp.K", "upper", "lower")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
  }
}

legend("bottomright", inset=.02, legend=c("Uniform", "Gamma", "Exponential"),
                bg="white",
                cex=0.8,
                pch=21,
                pt.bg=c("darkcyan", "grey", "orange"),
                lwd=1.5,
                lty=1,
                col=c("darkcyan", "grey", "orange"))

dev.off()

####################
#error - MSE

colors <- c("grey", "orange", "darkcyan", "lightblue", "#4D4D4D")
strategies <- c("gamma.st.05.01", "exponential.st.5", "uniform", "known.M", "full")

error <- aggregate(MSE ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
error <- cbind(error[,1:4], error$MSE)
colnames(error)[-c(1:4)] <- c("MSE", "lower", "upper")

pdf(file = "fig_3_5_MSE.pdf", width = 8, height = 4)

par(mfrow = c(1, 2))
par(mai=c(0.8, 0.8, 0.5, 0.2))

counter <- 0

for (v in c(0.1, 0.5)){
  for (A in 1){
    counter <- counter + 1
    
    plot(c(), c(), 
         ylim=c(0,0.4), 
         xlim=c(1, 8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1.1, 
         font.main = 1, 
         frame.plot=T, 
         xlab="", 
         ylab = "Mean Squared Error",
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    for (h in (1:4)/10){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
    }
    
    for (i in 1:length(strategies)){
      data.subset <- error[error$strategy==strategies[i] & error$alpha==A & error$var.exp==v, c("alleles", "MSE", "upper", "lower")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
  }
}

legend("topleft", inset=.02, legend=c("Uniform", "Gamma", "Exponential", "Full", "Known"),
                bg="white",
                cex=0.8,
                pch=21,
                pt.bg=colors[c(1:3,5,4)],
                lwd=1.5,
                lty=1,
                col=colors[c(1:3,5,4)])

dev.off()

####################
#accuracy - (0-1) - exponential

colors <- c("orange", "blue", "red", "purple")
strategies <- c("exponential.st.5", "exponential.st.5.tree", "exponential.st.5.tree.wrong", "exponential.st.5.tree.miss")
accuracy <- aggregate(accuracy ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=mean)

for (A in 1){
  pdf(file = paste0("fig_3_1_e_", A, "_accuracy.pdf"), width = 8, height = 4)
  
  par(mfrow = c(1, 2))
  par(mai=c(0.8, 0.8, 0.5, 0.2))
  
  counter <- 0
  
  for (v in c(0.1,0.5)){
    counter <- counter + 1
    
    plot(c(), c(), 
         ylim=c(0,1), 
         xlim=c(1,8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1.1, 
         font.main = 1, 
         frame.plot=T,
         xlab="", 
         ylab = "0-1 Accuracy",
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    for (h in (1:5)/5){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
    }
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "accuracy")]
      data.subset <- cbind(data.subset, t(sapply(data.subset[,2], binomial.prop.ci)))
      colnames(data.subset)[3:4] <- c("upper", "lower")
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
    
    if (v==0.1){
      legend("topright", inset=.02, legend=c("CRP", "Tree", "Misspecified", "Incorrect"),
             bg="white",
             cex=0.8,
             pch=21,
             pt.bg=colors,
             lwd=1.5,
             lty=1,
             col=colors)
    }
  }
}

dev.off()

####################
#accuracy - posterior.M - exponential

accuracy <- aggregate(posterior.M ~ strategy + alleles + alpha + var.exp, data=results.full[results.full$strategy %in% strategies,], FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
accuracy <- cbind(accuracy[,1:4], accuracy$posterior.M)
colnames(accuracy)[-c(1:4)] <- c("posterior.M", "lower", "upper")

for (A in 1){
  pdf(file = paste0("fig_3_2_e_", A, "_posterior_M.pdf"), width = 8, height = 4)
  
  par(mfrow = c(1, 2))
  par(mai=c(0.8, 0.8, 0.5, 0.2))
  
  counter <- 0
  
  for (v in c(0.1, 0.5)){
    counter <- counter + 1
    
    plot(c(), c(), 
         ylim=c(0,1), 
         xlim=c(1,8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1, 
         font.main = 1, 
         frame.plot=T,
         xlab="", 
         ylab = "Posterior Certainty",
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    for (h in (1:5)/5){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
    }
    
    for (i in 1:length(strategies)){
      data.subset <- accuracy[accuracy$strategy==strategies[i] & accuracy$alpha==A & accuracy$var.exp==v, c("alleles", "posterior.M", "lower", "upper")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
  }
  
  legend("topright", inset=.02, legend=c("CRP", "Tree", "Misspecified", "Incorrect"),
                  bg="white",
                  cex=0.8,
                  pch=21,
                  pt.bg=colors,
                  lwd=1.5,
                  lty=1,
                  col=colors)
  
  dev.off()
}

####################
#error - MSE - exponential

colors <- c("orange", "blue", "red", "purple", "lightblue", "#4D4D4D")
strategies <- c("exponential.st.5", "exponential.st.5.tree", "exponential.st.5.tree.wrong", "exponential.st.5.tree.miss", "known.M", "full")

error <- aggregate(MSE ~ strategy + alleles + alpha + var.exp, data=results.full, FUN=function(x){unlist(t.test(x)[c("estimate", "conf.int")])})  
error <- cbind(error[,1:4], error$MSE)
colnames(error)[-c(1:4)] <- c("MSE", "lower", "upper")

counter=0

for (A in 1){
  pdf(file = paste0("fig_3_5_e_", A, "_MSE.pdf"), width = 8, height = 4)
  
  par(mfrow = c(1, 2))
  par(mai=c(0.8, 0.8, 0.5, 0.2))
  
  for (v in c(0.1, 0.5)){
    counter <- counter+1
    plot(c(), c(), 
         ylim=c(0,0.4), 
         xlim=c(1, 8), 
         las=1, 
         cex.lab=1, 
         cex.axis=1, 
         cex.main=1.1, 
         font.main = 1, 
         frame.plot=T,
         xlab="", ylab = "Mean Squared Error",
         main=paste0("QTL Effect Size: ", v))
    
    title(xlab="Number of Alleles", line=2.5)
    
    for (h in (1:4)/10){
      abline(h=h, lty=2, col=scales::alpha("black", 0.3), lwd=0.7)
    }
    
    for (i in 1:length(strategies)){
      data.subset <- error[error$strategy==strategies[i] & error$alpha==A & error$var.exp==v, c("alleles", "MSE", "upper", "lower")]
      
      polygon(c(data.subset$alleles,rev(data.subset$alleles)), c(data.subset$lower,rev(data.subset$upper)), col=scales::alpha(colors[i], alpha=0.3), density=NA)
      lines(data.subset$alleles, data.subset[,2], lwd=1.5, col=colors[i])
      points(data.subset$alleles, data.subset[,2], col=colors[i], pch=19)
    }
    
    put.fig.letter(LETTERS[counter], "topleft", font=2, cex=1.2, offset=c(0.005, -0.005))
  }
}

legend("topleft", inset=.02, legend=c("CRP", "Tree", "Misspecified", "Incorrect", "Full", "Known"),
                bg="white",
                cex=0.8,
                pch=21,
                pt.bg=colors[c(1:4,6,5)],
                lwd=1.5,
                lty=1,
                col=colors[c(1:4,6,5)])

dev.off()
