##---------------------------------------- ##
## Parametric 2PL - constraints on abilities ----
##----------------------------------------##
code <- nimbleCode({
  for(i in 1:I) {
    for(j in 1:N) {
      y[j, i] ~ dbern(pi[j, i])

      pi[j,i] <- delta[i] + (1 - delta[i]) * linearReg[j, i]
      logit(linearReg[j, i]) <- lambda[i]*(eta[j] - beta[i])

    }
  }  
  
  for(i in 1:I) {
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    beta[i] ~ dnorm(0,  var = 3)
    delta[i] ~ dbeta(4, 12)
  } 
  
  for(j in 1:N) {
    eta[j] ~ dnorm(0, 1)
  }

  # ## dummy nodes to track log porbability and log likelihood
  myLogProbAll   ~ dnorm(0,1)
  myLogProbSome  ~ dnorm(0,1)
  myLogLik       ~ dnorm(0,1)

})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(beta   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1))

inits$lambda <- exp(inits$log_lambda)

monitors <- c("beta", "lambda", "delta")
