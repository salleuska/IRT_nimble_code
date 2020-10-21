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

code2PL <- nimbleCode({

  for(j in 1:N) {
    y[j, 1:I] ~ dBernoulliVector(prob = pi[j, 1:I])
    logit(pi[j, 1:I]) <-  lambda[1:I]*(eta[j] - beta[1:I])
  }

  for(i in 1:I) {
    beta.tmp[i] ~ dnorm(0, var = 3)  
    logLambda.tmp[i] ~ dnorm(0.5, var = 0.5)   
  }
  
  log(lambda[1:I]) <- logLambda.tmp[1:I] - mean(logLambda.tmp[1:I])
  beta[1:I] <- beta.tmp[1:I] - mean(beta.tmp[1:I])

  for(j in 1:N) {
    eta[j] ~ dnorm(mu, sd = sd.eta)
  }  

  mu ~ dnorm(0, var = 3)
  sd.eta ~ dunif(0.0001, 10)

})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(beta.tmp   = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1), 
              sd.eta = 1, mu = 0)


monitors <- c("beta", "lambda", "sd.eta", "mu")
