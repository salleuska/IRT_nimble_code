##----------------------------------------##
code2PL <- nimbleCode({
  for(i in 1:I) {
    for(j in 1:N) {
      y[j, i] ~ dbern(pi[j, i])
      logit(pi[j, i]) <-  lambda[i]*(eta[j] - beta[i])
    }
  }  
  
  for(i in 1:I) { 
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    beta[i] ~ dnorm(0,  var = 3)
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


constants <- list(I= dim(data$y)[2], N = dim(data$y)[1], M = 50)

inits <- list(beta       = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              nu1 = 2.01, nu2 = 1.01, s2_mu = 2, 
              alpha   = 1, a = 1, b = 3)

inits$lambda <- exp(inits$log_lambda)

monitors <- c("beta", "lambda",  "zi", "muTilde", "s2Tilde", "alpha")
