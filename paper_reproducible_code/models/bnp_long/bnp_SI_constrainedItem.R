##---------------------------------------- ##
code2PL <- nimbleCode({
  
  for(i in 1:NTot) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <-  lambda[item[i]]*eta[student[i]] + gamma[item[i]]
  }
  
  for(i in 1:I) {
    gamma.tmp[i] ~ dnorm(0, var = 3)  
    logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)   
  }
  
  m.logLambda <- sum(logLambda.tmp[1:I])/I
  m.gamma <- sum(gamma.tmp[1:I])/I
  
  for(i in 1:I) {
    gamma[i] <- gamma.tmp[i] - m.gamma
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
})


constants <- list(NTot= length(data$y),
                  I = length(unique(alldata$item)), 
                  N = length(unique(alldata$id)), 
                  student = alldata$id, 
                  item = alldata$item, 
                  M = 50)

inits <- list(gamma.tmp       = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1),
              nu1 = 2.01, nu2 = 1.01, s2_mu = 2, 
              alpha   = 1, a = 1, b = 3)

inits$lambda <- exp(inits$logLambda.tmp - mean(inits$logLambda.tmp))
inits$gamma <- inits$gamma.tmp - mean(inits$gamma.tmp)


monitors <- c("gamma", "lambda", "eta", "zi", "muTilde", "s2Tilde", "alpha")

