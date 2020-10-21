##------------------------------------------------------##
## Function to use  different samplers in compareMCMCs
##------------------------------------------------------##
source("util/customSamplers.R")
## centered sampler only for gamma parametrization 

configure_centeredSampler <- function(model) {
  conf <- configureMCMC(model)

  conf$removeSamplers("log_lambda")
  
  conf$addSampler(type = 'centered',
       target = c('log_lambda', 'gamma'),
       control = list(nodesToCenter = 'eta', scale = 0.1, adaptive = TRUE))
  conf
}

## Paired sampler

assignPairedSampling <- function(mcmcConf, var, scale = 0.5, removeExistingSamplers = FALSE) {
    nodeNames <- mcmcConf$model$expandNodeNames(var)
    
    if(removeExistingSamplers)    mcmcConf$removeSamplers(var)
    nNodes <- length(nodeNames)
    for(i in 1:(nNodes-1)) {
       mcmcConf$addSampler(type = "RWpaired", target = nodeNames[i:(i+1)], control = list(scale = scale, adaptive = TRUE))
    }
}


configure_pairedSamplerBeta <- function(model) {
  conf <- configureMCMC(model)

  conf$removeSamplers("log_lambda")
  assignPairedSampling(conf, "log_lambda", scale = 0.1)
  conf$removeSamplers("beta")
  assignPairedSampling(conf, "beta", scale = 1)
  
  conf
}


configure_pairedSamplerGamma <- function(model) {
  conf <- configureMCMC(model)

  conf$removeSamplers("log_lambda")
  assignPairedSampling(conf, "log_lambda", scale = 0.1)
  conf$removeSamplers("gamma")
  assignPairedSampling(conf, "gamma", scale = 1)
  
  conf
}