################################################################
#########  FUNCTIONS
################################################################
## posterior rescaling
## account for eta thinned

posteriorRescalingGamma <- function(samples, samples2 = NULL, thinEta = 1,  rescale = TRUE){
  
  lambdaSamp <- samples[, grep("^lambda", colnames(samples))]
  gammaSamp <- samples[, grep("gamma", colnames(samples))]
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
  betaSamp <- samples[, grep("beta", colnames(samples))]
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


# plotItemsEstimates <- function(samples){
# if("betaSamp" %in% names(samples)){
#   par(mfrow = c(1, 3))
#   plot(apply(samples$lambdaSamp, 2, mean), lambda0.disc, 
#        xlab =  TeX('Estimated $\\hat{\\lambda}$'),
#        ylab = TeX('True $\\lambda$'), 
#        main = "Discrimination parameters")
#   abline(c(0,1))
#   plot(apply(- samples$beta*samples$lambdaSamp,2, mean), -beta0.diff*lambda0.disc, 
#        xlab =  TeX('Estimated $\\hat{\\gamma}$'),
#        ylab = TeX('True $\\gamma$'), 
#        main = "Difficulty parameters - reparametrization")
#   abline(c(0,1))
#   plot(apply(samples$beta, 2, mean), beta0.diff,
#        xlab =  TeX('Estimated $\\hat{\\beta} = -\\hat{\\gamma}/\\hat{\\lambda}$'),
#        ylab = TeX('True $\\beta$'),
#        main = "Difficulty parameters")
#   abline(c(0,1))
# } else {
#   par(mfrow = c(1, 3))
#   plot(apply(samples$lambdaSamp, 2, mean), lambda0.disc, 
#        xlab =  TeX('Estimated $\\hat{\\lambda}$'),
#        ylab = TeX('True $\\lambda$'), 
#        main = "Discrimination parameters")
#   abline(c(0,1))
#   plot(apply(samples$gammaSamp, 2, mean), -beta0.diff*lambda0.disc, 
#        xlab =  TeX('Estimated $\\hat{\\gamma}$'),
#        ylab = TeX('True $\\gamma$'), 
#        main = "Difficulty parameters - reparametrization")
#   abline(c(0,1))
#   plot(apply(-samples$gammaSamp/samples$lambdaSamp, 2, mean), beta0.diff,
#        xlab =  TeX('Estimated $\\hat{\\beta} = -\\hat{\\gamma}/\\hat{\\lambda}$'),
#        ylab = TeX('True $\\beta$'),
#        main = "Difficulty parameters")
#   abline(c(0,1))  
#  }
# }


# makeTraceplotReport <- function(samples, filename, title_Traceplots){
# rmarkdown::render("Traceplots.Rmd", 
#                   output_file = paste0(dirOutput, "Traceplot_", filename , ".html"), 
#                   params = list(nimbleSamples = samples,
#                                 set_title = paste0(title_Traceplots, "-", filename)))  
# }

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