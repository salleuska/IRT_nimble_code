##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
## posterior rescaling
## account for abilities parameters thinned

posteriorRescalingGamma <- function(samples, samples2 = NULL, thinEta = 1,  rescale = TRUE){
  
  lambdaSamp <- samples[, grep("^lambda", colnames(samples))]
  gammaSamp <- samples[, grep("gamma\\[", colnames(samples))]
  otherParSamp <- samples[, -grep("(lambda|gamma|^eta)", colnames(samples))]

  if(is.null(samples2)){
    etaSamp <- samples[, grep("^eta", colnames(samples))]
  } else {
    etaSamp <- samples2[, grep("^eta", colnames(samples2))]   
  }

  scaleShiftEta <- NULL
  locationShiftEta <- NULL

  if(rescale){
    
    nItems <- dim(lambdaSamp)[2]
    nSamp <- dim(lambdaSamp)[1]
    
    locationSamp <- sapply(1:nSamp, function(x) sum(gammaSamp[x,])/sum(lambdaSamp[x,]))
    scaleSamp <- apply(lambdaSamp, 1, function(x) prod(x)^(-1/nItems))
    
    gammaSamp <- t(sapply(1:nSamp, function(x) gammaSamp[x, ] - lambdaSamp[x,]*locationSamp[x]))
    lambdaSamp <- t(apply(lambdaSamp, 1, function(x) x*(prod(x)^(-1/nItems)) ))
    
    indicesEta <- seq(thinEta, nSamp, by = thinEta)

    etaSamp2 <- matrix(0, NROW(etaSamp), NCOL(etaSamp))
    for(i in 1:dim(etaSamp)[1]){
      etaSamp2[i, ] <- (etaSamp[i, ] + locationSamp[indicesEta[i] ])/scaleSamp[indicesEta[i]]
    }
    colnames(etaSamp2) <- colnames(etaSamp) 
    
    etaSamp <- etaSamp2

    scaleShiftEta = scaleSamp[indicesEta]
    locationShiftEta = locationSamp[indicesEta]
    
  }
  
  out <- list(lambdaSamp = lambdaSamp, 
              gammaSamp = gammaSamp, 
              etaSamp = etaSamp, 
              scaleShiftEta = scaleShiftEta,
              locationShiftEta = locationShiftEta,
              otherParSamp = otherParSamp)
  out 
}



posteriorRescalingBeta <- function(samples, samples2 = NULL, thinEta = 1, rescale = TRUE){
  
  lambdaSamp <-samples[, grep("^lambda", colnames(samples))]
  betaSamp <- samples[, grep("beta\\[", colnames(samples))]
  otherParSamp <- samples[, -grep("(lambda|beta|^eta)", colnames(samples))]

  if(is.null(samples2)){
    etaSamp <- samples[, grep("^eta", colnames(samples))]
  } else {
    etaSamp <- samples2[, grep("^eta", colnames(samples2))]   
  }
  
  scaleShiftEta <- NULL
  locationShiftEta <- NULL

  if(rescale){
    
    nItems <- dim(lambdaSamp)[2]
    nSamp <- dim(lambdaSamp)[1]
    
    locationSamp <- sapply(1:nSamp, function(x) sum(betaSamp[x,])/nItems)
    scaleSamp <- apply(lambdaSamp, 1, function(x) prod(x)^(-1/nItems))
    
    betaSamp <- t(sapply(1:nSamp, function(x) (betaSamp[x, ] - locationSamp[x])/scaleSamp[x] ))
    lambdaSamp <- t(apply(lambdaSamp, 1, function(x) x*(prod(x)^(-1/nItems)) ))
    
    indicesEta <- seq(thinEta, nSamp, by = thinEta)

    etaSamp2 <- matrix(0, NROW(etaSamp), NCOL(etaSamp))
    for(i in 1:dim(etaSamp)[1]){
      etaSamp2[i, ] <- (etaSamp[i, ] - locationSamp[indicesEta[i]])/scaleSamp[indicesEta[i]]
    }
    colnames(etaSamp2) <- colnames(etaSamp) 
    etaSamp <- etaSamp2

    scaleShiftEta = scaleSamp[indicesEta]
    locationShiftEta = locationSamp[indicesEta]

  }

  
  out <- list(lambdaSamp = lambdaSamp, 
              betaSamp = betaSamp, 
              etaSamp = etaSamp, 
              scaleShiftEta = scaleShiftEta,
              locationShiftEta = locationShiftEta,
              otherParSamp = otherParSamp)
  out 
}


getNewCompareMCMCObjBeta <- function(resCompareMCMC, rescaledSamples, sampleName) {
  resCompareMCMC[[1]]$samples <-  cbind(rescaledSamples$lambdaSamp,rescaledSamples$betaSamp)
  
  resCompareMCMC[[1]]$MCMC <- sampleName
  
  clearMetrics(resCompareMCMC)
  addMetrics(resCompareMCMC) # use default metrics
  combineMetrics(resCompareMCMC)
  
  resCompareMCMC
}

getNewCompareMCMCObjGamma <- function(resCompareMCMC, rescaledSamples, sampleName) {
  resCompareMCMC[[1]]$samples <-  cbind(rescaledSamples$lambdaSamp,rescaledSamples$gammaSamp)
  
  resCompareMCMC[[1]]$MCMC <- sampleName
  
  clearMetrics(resCompareMCMC)
  addMetrics(resCompareMCMC) # use default metrics
  combineMetrics(resCompareMCMC)
  
  resCompareMCMC
}