  ##---------------------------------------- ##
## Parametric 2PL - unconstrained ----
##----------------------------------------##
code2PL <- nimbleCode({
  for(i in 1:NTot) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <-  lambda[item[i]]*(eta[student[i]] - beta[item[i]])
  }
  
  for(i in 1:I) {
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    beta[i] ~ dnorm(0,  var = 3)
  } 
  
  
  for(j in 1:N) {
    eta[j] ~ dnorm(mu, var = s2.eta)
  }  

  mu ~ dnorm(0, var = 3)
  s2.eta ~ dinvgamma(2.01, 1.01)
})

constants <- list(NTot= length(data$y)[1],
                  I = length(unique(alldata$item)), 
                  N = length(unique(alldata$id)), 
                  student = alldata$id, 
                  item = alldata$item)

inits <- list(beta   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              s2.eta = 1, mu = 0)

inits$lambda <- exp(inits$log_lambda)

monitors <- c("beta", "lambda", "eta", "s2.eta", "mu")
