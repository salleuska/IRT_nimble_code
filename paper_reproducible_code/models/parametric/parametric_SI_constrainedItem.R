##---------------------------------------- ##
## Parametric 2PL - constrained ----
##----------------------------------------##

dBernoulliVector <- nimbleFunction(
  run = function(x    = double(1), 
                 prob = double(1), 
                 log  = integer(0, default = 0)) {

    returnType(double(0))
    logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  }
)

rBernoulliVector <- nimbleFunction(
  run = function(n    = integer(0), 
                 prob = double(1)) {

    returnType(double(1))
    n = length(prob)
    return(rbinom(n, size = 1, prob = prob))
  }
)

code <- nimbleCode({
 
  for(j in 1:N) {
    y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])
    logit(pi[j, 1:I]) <- lambda[1:I]*eta[j] + gamma[1:I]
  }

  for(i in 1:I) {
    gamma.tmp[i] ~ dnorm(0, var = 3)  
    logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)   
  }
  
  log(lambda[1:I]) <- logLambda.tmp[1:I] - mean(logLambda.tmp[1:I])
  gamma[1:I] <- gamma.tmp[1:I] - mean(gamma.tmp[1:I])

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

inits <- list(gamma.tmp   = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1),
              s2.eta = 1, mu = 0)


monitors <- c("gamma", "lambda", "s2.eta", "mu")
