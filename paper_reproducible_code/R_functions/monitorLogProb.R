##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
#############
## Custom sampler to monitor posterior logProbabilities
############# 

## sampler to monitor posterior log-probabilities 
logProb_summer <- nimbleFunction(
	name = 'logProb_summer',
	contains = sampler_BASE,
    setup = function(model, mvSaved, target,  control) {
    	if(!is.null(control$nodeList)){
    		nodes <- model$expandNodeNames(control$nodeList)
    	} else {
        ## use all nodes in the model if the user does not provide
        ## a list of nodes to be monitored
  			nodes <- model$getNodeNames()
    	}
      nodes <- nodes[nodes != target]

    },
    run = function() {
        model[[target]] <<- model$getLogProb(nodes)  
        
    copy(from = model, to = mvSaved, 
    	row = 1, nodes = target, logProb = TRUE)
   }, 
   methods = list( reset = function () {} )
)
 
# logLik_summer <- nimbleFunction(
# 	name = 'logLik_summer',
# 	contains = sampler_BASE,
#     setup = function(model, mvSaved, target, control) {
#           dataNodes <- model$getNodeNames(dataOnly = TRUE)
#           dataNodes <- dataNodes[dataNodes != target]
#     },
#     run = function() {
#         model[[target]] <<- model$getLogProb(dataNodes)  

#         copy(from = model, to = mvSaved, 
#     	row = 1, nodes = target, logProb = TRUE)

#    }, 
#    methods = list( reset = function () {})
# )
 

# # monitor logProb_sum
# conf$removeSampler("logProb_sum")
# conf$addSampler("logProb_sum", type =  "logProb_summer", 
# 				control = list(nodeList = c("beta", "lambda", "eta")) )

# data <- list(y = readRDS("data/simulation_bimodal.rds"))

# ###------------------------------------------------ ##
# ## Parametric 2PL - constraints on abilities ----
# ##------------------------------------------------##
# code2PL <- nimbleCode({
#   for(i in 1:I) {
#     for(j in 1:N) {
#       y[j, i] ~ dbern(pi[j, i])
#       logit(pi[j, i]) <-  lambda[i]*(eta[j] - beta[i])
#     }
#   }  
  
#   for(i in 1:I) {
#     log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
#     beta[i] ~ dnorm(0,  var = 3)
#   } 
  
#   for(j in 1:N) {
#     eta[j] ~ dnorm(0, 1)
#   }

#   # ## dummy nodes to track log porbability and log likelihood
#   myLogProbAll  ~ dnorm(0,1)
#   myLogProbSome ~ dnorm(0,1)
#   myLogLik      ~ dnorm(0,1)
# })

# constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

# inits <- list(beta   = rnorm(constants$I, 0, 1),
#               log_lambda = runif(constants$I, -1, 1))

# inits$lambda <- exp(inits$log_lambda)

# monitors <- c("beta", "lambda")


# model <- nimbleModel(code  = code2PL,
# 					data 			= data,  
# 					constants	= constants,
# 					inits 			= inits, 
# 					calculate 	= FALSE)

# monitors <- c(monitors, "myLogProbAll", "myLogProbSome", "myLogLik")
# # monitors <- c(monitors, "myLogProbAll")

# mcmcConf <- configureMCMC(model, monitors = monitors)

# mcmcConf$removeSampler("myLogLik")
# mcmcConf$addSampler("myLogLik", type =  "logProb_summer", 
#   control = list(nodeList = c("y")))

# mcmcConf$removeSampler("myLogProbAll")
# mcmcConf$addSampler("myLogProbAll", type =  "logProb_summer")

# mcmcConf$removeSampler("myLogProbSome")
# mcmcConf$addSampler("myLogProbSome", type =  "logProb_summer", 
#   control = list(nodeList = c("beta", "lambda", "eta")))
# mcmcConf

# mcmc <- buildMCMC(mcmcConf)	

# Cmodel <- compileNimble(model)
# Cmcmc <- compileNimble(mcmc, project = model)


# res <- runMCMC(Cmcmc, 
# 		   niter 	= 1000,  
# 		   nburnin  = 100, 
# 		   setSeed  = 32)

# colnames(res)
# res[, "myLogProbAll"]
# res[, "myLogProbSome"]
# res[, "myLogLik"]













