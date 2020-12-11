##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
library(nimble)
source("R_functions/rescalingFunctions.R")
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --resFileName
## --outDir
## resFileName="output/posterior_samples/simulation_bimodal/parametric/parametric_SI_unconstrained.rds"
## args <- list(resFileName, outDir)
##-----------------------------------------#
if(is.null(args$outDir)) outDir <- "output/posterior_samples_elaborated/" else dir <- args$outDir

data <- strsplit(resFileName, "\\/|.rds")[[1]][3]

fileName <- strsplit(resFileName, "\\/|.rds")[[1]][5]

## parameterization
param <- strsplit(basename(fileName), "\\_|.rds")[[1]][2]
## model type
modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]
## constraint
constraint <- strsplit(basename(fileName), "\\_|.rds")[[1]][3]

## read objects
resObj <- readRDS(resFileName)

##-------------------------------------------------------##
## rescale posterior samples to common parameterization
## IRT constrainedItem
##-------------------------------------------------------##

## set flag to true by default
if(constraint == "constrainedItem") rescale <- FALSE else rescale <- TRUE 

if(constraint == "stan") { 
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

outDirResults <- paste0(outDir, data, "/", modelType)
dir.create(file.path(outDirResults), recursive = TRUE, showWarnings = FALSE)

saveRDS(modelRes, file = paste0(outDirResults, "/" , fileName, ".rds"))

##-----------------------------------------#
## Get efficency results 
##-----------------------------------------#
ess_coda   <- NA
ess_multi  <- NA

xx <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
            modelRes[grepl("beta", names(modelRes))][[1]])

## compute ess for item parameters using different packages
ess_coda <- min(coda::effectiveSize(xx))
ess_multi <- mcmcse::multiESS(xx) ## NaN negative eigenvalue

compilationTime <- resObj$compilationTime[3]
runningTime <- resObj$runningTime[3]
samplingTime <- 0

## if parametric save also sampling time
if(modelType == "parametric"){ 
  if(constraint == "stan") { 
    samplingTime <- resObj$samplingTime
  } else {
    percBurnin <- resObj$MCMCcontrol$nburnin/resObj$MCMCcontrol$niter
    samplingTime <- runningTime * (1 - percBurnin)
  }

}

outDirTime <- paste0("output/mcmc_time/", data)
dir.create(file.path(outDirTime), recursive = TRUE, showWarnings = FALSE)

outFile <- paste0(outDirTime, "/", modelType, "efficiency.txt")
row <- cbind(fileName, ess_coda, ess_multi, compilationTime, runningTime, samplingTime)

if(!file.exists(outFile)){
	cat(colnames(row), "\n", file = outFile)
}
# append row
cat(row, "\n", file = outFile, append = TRUE)
# ############################################################################