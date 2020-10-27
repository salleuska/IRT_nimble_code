
##---------------------------------------- ##
## Parametric 2PL - unconstrained ----
##----------------------------------------##
code2PL <- nimbleCode({
  for(i in 1:I) {
    for(j in 1:N) {
      y[j, i] ~ dbern(pi[j, i])
      logit(pi[j, i]) <-  lambda[i]*eta[j] + gamma[i]
    }
  }  
  
  for(i in 1:I) { 
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    gamma[i] ~ dnorm(0,  var = 3)
  } 
  
  
  for(j in 1:N) {
    eta[j] ~ dnorm(mu, sd = sd.eta)
  }  

  mu ~ dnorm(0, var = 3)
  sd.eta ~ dunif(0.0001, 10)

})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(gamma   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              sd.eta = 1, mu = 0)

inits$lambda <- exp(inits$log_lambda)

monitors <- c("gamma", "lambda", "sd.eta", "mu")
