
##---------------------------------------- ##
## Parametric 2PL - unconstrained ----
##----------------------------------------##
code <- nimbleCode({
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
    eta[j] ~ dnorm(mu, var = s2.eta)
  }  

  mu ~ dnorm(0, var = 3)
  s2.eta ~ dinvgamma(2.01, 1.01)


  ## dummy nodes to track log porbability and log likelihood
  myLogProbAll   ~ dnorm(0,1)
  myLogProbSome  ~ dnorm(0,1)
  myLogLik       ~ dnorm(0,1)

})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(gamma   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              s2.eta = 1, mu = 0)

inits$lambda <- exp(inits$log_lambda)

monitors <- c("gamma", "lambda", "s2.eta", "mu")
