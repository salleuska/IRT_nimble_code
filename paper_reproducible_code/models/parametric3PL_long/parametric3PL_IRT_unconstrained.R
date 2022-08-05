##---------------------------------------- ##
## 3PL model - parametric
##----------------------------------------##
## lambda[i] - discrimination parameter
## beta[i]   - difficulty paramter
## delta[i]  - guessing parameter
## theta[j]  - individial ability

code <- nimbleCode({
  for(i in 1:NTot) {
      y[i] ~ dbern(pi[i])
      pi[i] <- delta[item[i]] + (1 - delta[item[i]]) * linearReg[i]
      logit(linearReg[i]) <-  lambda[item[i]]*(eta[student[i]] - beta[item[i]])
  }
    
   
  for(i in 1:I) {
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    beta[i] ~ dnorm(0,  var = 3)
     delta[i] ~ dbeta(4, 12)
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

inits <- list(beta   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1),
              delta  = rbeta(constants$I, 4, 12),
              s2.eta = 1, mu = 0)

inits$lambda <- exp(inits$log_lambda)

monitors <- c("beta", "lambda", "delta", "s2.eta", "mu")
