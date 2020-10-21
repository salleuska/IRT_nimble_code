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
    eta[j] ~ dnorm(mu, sd = sd.eta)
  }  

  mu ~ dnorm(0, var = 3)
  sd.eta ~ dunif(0.0001, 10)
})

constants <- list(NTot= length(data$y)[1],
                  I = length(unique(alldata$item)), 
                  N = length(unique(alldata$id)), 
                  student = alldata$id, 
                  item = alldata$item)

inits <- list(beta   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              sd.eta = 1, mu = 0)

inits$lambda <- exp(inits$log_lambda)

monitors <- c("beta", "lambda", "eta", "sd.eta", "mu")
