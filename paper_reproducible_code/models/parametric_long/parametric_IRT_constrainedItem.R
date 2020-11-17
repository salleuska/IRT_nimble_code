##---------------------------------------- ##
## Parametric 2PL - constrained ----
##----------------------------------------##
code2PL <- nimbleCode({
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

inits <- list(beta.tmp   = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1), 
              sd.eta = 1, mu = 0)

monitors <- c("beta", "lambda", "eta", "sd.eta", "mu")