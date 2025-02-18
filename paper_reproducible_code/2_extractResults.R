##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
## This script post-process the posterior samples
## and calculate different metrics (ess, WAIC)
##-----------------------------------------#
library(nimble)
source("R_functions/rescalingFunctions.R")
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --outDir       directory for output
## --resFileName  file name containing results
##-----------------------------------------#
if(is.null(args$outDir)) outDir <- "output/posterior_samples_elaborated/" else dir <- args$outDir

listLength <- length(strsplit(args$resFileName, "\\/|.rds")[[1]])

data <- strsplit(args$resFileName, "\\/|.rds")[[1]][listLength -2]

fileName <- strsplit(args$resFileName, "\\/|.rds")[[1]][listLength]

## parameterization
param <- strsplit(basename(fileName), "\\_|.rds")[[1]][2]
## model type
modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]
## constraint
constraint <- strsplit(basename(fileName), "\\_|.rds")[[1]][3]

## read objects
resObj <- readRDS(args$resFileName)


if(grepl("stan", constraint)) { 
  thinEta <-1 
  nSamp <- resObj$MCMCcontrol$niter - resObj$MCMCcontrol$nwarmup
  indicesEta <- seq(from = thinEta,  
            to = nSamp, 
            by = thinEta)
} else {
  ## Check and store thinning info for eta
  thinEta <- resObj$MCMCcontrol$thin2
  nSamp <- resObj$MCMCcontrol$niter - resObj$MCMCcontrol$nburnin

  indicesEta <- seq(from = thinEta,  
            to = nSamp, 
            by = thinEta)  
}


##-------------------------------------------------------##
## rescale posterior samples to common parameterization
## IRT constrainedItem
##-------------------------------------------------------##

## set flag to true by default
if(constraint == "constrainedItem") rescale <- FALSE else rescale <- TRUE 

if(grepl("stan", constraint)) { 
  modelRes <- posteriorRescalingBeta(samples  = resObj$samples[, -grep("^eta", colnames(resObj$samples))],
                                     samples2 = resObj$samples[, grep("^eta", colnames(resObj$samples))],
                                     thinEta  = 1, 
                                     rescale  = rescale)
} else if(param == "IRT" ){
  modelRes <- posteriorRescalingBeta(samples  = resObj$samples$samples,
                                     samples2 = resObj$samples$samples2,
                                     thinEta  = resObj$MCMCcontrol$thin2, 
                                     rescale  = rescale)
} else if(param == "SI") {
  modelRes <- posteriorRescalingGamma(samples  = resObj$samples$samples,
                                      samples2 = resObj$samples$samples2,
                                      thinEta  = resObj$MCMCcontrol$thin2, 
                                      rescale  = rescale)
  
  modelRes$betaSamp <- -modelRes$gammaSamp/modelRes$lambdaSamp
  modelRes$betaSamp <-  modelRes$betaSamp - apply(modelRes$betaSamp, 1, mean) 
  colnames(modelRes$betaSamp) <- gsub("gamma", "beta", colnames(modelRes$betaSamp))
}


## matrices for ESS evaluations
onlyItems <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
              modelRes[grepl("beta", names(modelRes))][[1]])

itemsAndAbility <- cbind(onlyItems[indicesEta, ], modelRes$etaSamp)


## for plots etc save thinned samples for abilities 
## and other related parameters

if(constraint != "stan" & thinEta == 1) {
  thinEta <- 10
  indicesEta <- seq(from = thinEta,  
            to = nSamp, 
            by = thinEta)  

  modelRes$etaSamp <- modelRes$etaSamp[indicesEta, ]
  modelRes$scaleShiftEta <- modelRes$scaleShiftEta[indicesEta]
  modelRes$locationShiftEta <- modelRes$locationShiftEta[indicesEta]
  modelRes$otherParSamp <- modelRes$otherParSamp[indicesEta, ]

}
outDirResults <- paste0(outDir, data, "/")
dir.create(file.path(outDirResults), recursive = TRUE, showWarnings = FALSE)

#######################
saveRDS(modelRes, file = paste0(outDirResults, "/" , fileName, ".rds"))
#######################
##-----------------------------------------#
## Get efficency results 
##-----------------------------------------#
## all parameters at same level (item + sampled abilities)

essCodaItems          <- NA
essCodaItemsAbility   <- NA

essCodaLogLik                <- NA
essCodaLogPostAll            <- NA
essCodaLogPostItemsAbility   <- NA

multiEssItemsAbility   <- NA


#######################
## compute mixing performance measures
essCodaItems <- min(coda::effectiveSize(onlyItems))

##compute and save ESS for parameters
essItemsAbility <- coda::effectiveSize(itemsAndAbility)

essCodaItemsAbility <- min(essItemsAbility)

if(!grepl("stan", constraint)) { 
  essCodaLogLik                <- coda::effectiveSize(modelRes$otherParSamp[, "myLogLik"])
  essCodaLogPostAll            <- coda::effectiveSize(modelRes$otherParSamp[, "myLogProbAll"])
  essCodaLogPostItemsAbility   <- coda::effectiveSize(modelRes$otherParSamp[, "myLogProbSome"])
}
#######################
itemsAndAbilityMultiESS <- itemsAndAbility[, !grepl("(beta\\[1\\])|(lambda\\[1\\])", colnames(itemsAndAbility))]
#######################
## Multivariate ESS using multivariate Batch mean(mBm) stimator 
## Vats et al Biometrika paper on multivariate ESS
## note: by defaul multiESS tries to get the lugsail estimator
## (i.e. starts from the mBm, applying some corrections)
## which does not work in this application

try(multiEssItemsAbility <- mcmcse::multiESS(itemsAndAbilityMultiESS, 
                                            method = "bm", r = 1, adjust = FALSE))


# xx <- mcmcse::multiESS(itemsAndAbilityMultiESS, 
#                                             method = "bm", r = 3)
# try(essItemsAbilityMCMCSE <- mcmcse::ess(itemsAndAbilityMultiESS, 
#                                             method = "bm", r = 1))

## Extract times

compilationTime <- resObj$compilationTime[3]
runningTime <- resObj$runningTime[3]
samplingTime <- 0

## TO plot distribtion of the ESS
#saveRDS(essItemsAbilityMCMCSE/runningTime, file = paste0(outDirResults, "/ESS_" , fileName, ".rds"))

#######################
## if parametric save also sampling time
if(modelType == "parametric"){ 
  if(grepl("stan", constraint)) { 
    samplingTime <- resObj$samplingTime
  } else {
    percBurnin <- resObj$MCMCcontrol$nburnin/resObj$MCMCcontrol$niter
    samplingTime <- runningTime * (1 - percBurnin)
  }

}

WAIC <- 0
#if(modelType == "parametric"){ 
if(modelType %in% c("parametric", "parametric3PL", "bnp3PL")){ 
  WAIC <- resObj$modelWAIC
} 
if(is.null(WAIC)){ WAIC <-0} 

outDirTime <- paste0("output/mcmc_time/", data)
dir.create(file.path(outDirTime), recursive = TRUE, showWarnings = FALSE)

outFile <- paste0(outDirTime, "/", modelType, "_efficiency.txt")


row <- cbind(fileName, 
              essCodaItems,
              essCodaItemsAbility,
              essCodaLogLik,
              essCodaLogPostAll,
              essCodaLogPostItemsAbility,
              multiEssItemsAbility,
              compilationTime,
              runningTime,
              samplingTime,
              WAIC)

if(!file.exists(outFile)){
	## if file does not exist create it and append row
  cat(colnames(row), "\n", file = outFile)
  # append row
  cat(row, "\n", file = outFile, append = TRUE)
} else {
  file <- read.table(outFile, header = T)

  if(row[, "fileName"] %in% file$fileName) {
    file[which(row[, "fileName"] == file$fileName), ] <- row[1, ]
  } else {    
    file <- rbind(file, row[1, ])
  }
  
  write.table(file, file = outFile, 
      col.names = TRUE, row.names = FALSE, quote = FALSE)

}

############################################################################

