dBernoulliVector <- nimbleFunction(
  run = function(x    = double(1), 
                 prob = double(1), 
                 log  = integer(0)) {

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
    
  ## CRP for clustering individual effects
  zi[1:N] ~ dCRP(alpha, size = N)
  alpha ~ dgamma(a, b)  
  ## Mixture component parameter drawn from the base measure
  for(j in 1:N) {
    eta[j] ~ dnorm(mu[j], var = s2[j])  
    mu[j] <- muTilde[zi[j]]                 
    s2[j] <- s2Tilde[zi[j]]   
  }

  for(m in 1:M) {
    muTilde[m] ~ dnorm(0, var = s2_mu)
    s2Tilde[m] ~ dinvgamma(nu1, nu2)
  }
})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1], M = 50)

inits <- list(beta.tmp       = rnorm(constants$I, 0, 1),
              logLambda.tmp  = runif(constants$I, -1, 1),
              nu1 = 2.01, nu2 = 1.01, s2_mu = 2, 
              alpha   = 1, a = 1, b = 3)

inits$lambda <- exp(inits$logLambda.tmp - mean(inits$logLambda.tmp))
inits$beta <- inits$beta.tmp - mean(inits$beta.tmp)

monitors <- c("beta", "lambda",  "zi", "muTilde", "s2Tilde", "alpha")