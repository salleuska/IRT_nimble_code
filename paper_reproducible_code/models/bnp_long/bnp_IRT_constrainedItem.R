##----------------------------------------##
code <- nimbleCode({
  
  for(i in 1:NTot) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <-  lambda[item[i]]*(eta[student[i]] - beta[item[i]])
  }

  for(i in 1:I) {
    beta.tmp[i] ~ dnorm(0, var = 3)  
    logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)   
  }
  
  m.logLambda <- sum(logLambda.tmp[1:I])/I
  m.beta <- sum(beta.tmp[1:I])/I
  
  for(i in 1:I) {
    beta[i] <- beta.tmp[i] - m.beta
    log(lambda[i]) <- logLambda.tmp[i] - m.logLambda
  } 
    
  ## CRP for clustering individual effects
  zi[1:N] ~ dCRP(alpha, size = N)
  alpha ~ dgamma(a, b)  
  ## Mixture component parameter drawn from the base measure
  for(j in 1:N) {
    eta[j] ~ dnorm(mu[j], var = s2[j])  
    mu[j] <- muTilde[zi[j]]                 
    s2[j] <- s2Tilde[zi[j]]   
  }

  for(m in 1:M) {
    muTilde[m] ~ dnorm(0, var = s2_mu)
    s2Tilde[m] ~ dinvgamma(nu1, nu2)
  }
  
  ## dummy nodes to track log probability and log likelihood
  myLogProbAll   ~ dnorm(0,1)
  myLogProbSome  ~ dnorm(0,1)
  myLogLik       ~ dnorm(0,1)

})

constants <- list(NTot= length(data$y),
                  I = length(unique(alldata$item)), 
                  N = length(unique(alldata$id)), 
                  student = alldata$id, 
                  item = alldata$item, 
                  M = 50)

inits <- list(beta.tmp       = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1),
              nu1 = 2.01, nu2 = 1.01, s2_mu = 2, 
              alpha   = 1, a = 1, b = 3)

inits$lambda <- exp(inits$logLambda.tmp - mean(inits$logLambda.tmp))
inits$beta <- inits$beta.tmp - mean(inits$beta.tmp)

monitors <- c("beta", "lambda", "zi", "muTilde", "s2Tilde", "alpha")



