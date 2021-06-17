##---------------------------------------- ##
## Parametric 2PL - constraints on abilities ----
##----------------------------------------##
code2PL <- nimbleCode({
  for(i in 1:NTot) {
    y[i] ~ dbern(pi[i])
    logit(pi[i]) <-  lambda[item[i]]*eta[student[i]] + gamma[item[i]]
  }
  
  for(i in 1:I) {
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    gamma[i] ~ dnorm(0,  var = 3)
  } 
  
  
  for(j in 1:N) {
    eta[j] ~ dnorm(0, 1)
  }  

})


constants <- list(NTot= length(data$y)[1],
                  I = length(unique(alldata$item)), 
                  N = length(unique(alldata$id)), 
                  student = alldata$id, 
                  item = alldata$item)

inits <- list(gamma   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1))

inits$lambda <- exp(inits$log_lambda)

monitors <- c("gamma", "lambda")
