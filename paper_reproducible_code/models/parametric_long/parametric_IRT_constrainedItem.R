##---------------------------------------- ##
## Parametric 2PL - constrained ----
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

constants <- list(NTot= length(data$y)[1],
                  I = length(unique(alldata$item)), 
                  N = length(unique(alldata$id)), 
                  student = alldata$id, 
                  item = alldata$item)

inits <- list(beta.tmp   = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1), 
              s2.eta = 1, mu = 0)

monitors <- c("beta", "lambda", "s2.eta", "mu")
