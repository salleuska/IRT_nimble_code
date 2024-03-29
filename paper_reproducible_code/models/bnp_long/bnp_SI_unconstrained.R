##---------------------------------------- ##
code <- nimbleCode({
  for(i in 1:NTot) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <-  lambda[item[i]]*eta[student[i]] + gamma[item[i]]
  }
  
  for(i in 1:I) { 
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    gamma[i] ~ dnorm(0,  var = 3)
  } 
  
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

  
  ## dummy nodes to track log porbability and log likelihood
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

inits <- list(gamma       = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              nu1 = 2.01, nu2 = 1.01, s2_mu = 2, 
              alpha  = 0.2, a = 1, b = 3)

inits$lambda <- exp(inits$log_lambda)


monitors <- c("gamma", "lambda",  "zi", "muTilde", "s2Tilde", "alpha")
