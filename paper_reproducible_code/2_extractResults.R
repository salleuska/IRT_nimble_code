################################################################
library(nimble)
source("R_functions/rescalingFunctions.R")
################################################################
#########  SETTINGS
################################################################
args <- R.utils::commandArgs(asValue=TRUE)
## --dirResults=
## --dirOutput=
###############
# resFileName="simulation_bimodal/bnp/res/2PL_unconstrained_bnp_beta.rds"
# dirOutput="simulation_bimodal/bnp/elaboratedRes/"
# args <- list(resFileName= resFileName, dirOutput = dirOutput)

# resFileName="data_health/bnp/res/2PL_unconstrained_bnp_beta.rds"
# dirOutput="data_health/bnp/elaboratedRes/"
# args <- list(resFileName= resFileName, dirOutput = dirOutput)

# resFileName="data_timss/bnp/res/2PL_unconstrained_bnp_beta.rds"
# dirOutput="data_timss/bnp/elaboratedRes/"
# args <- list(resFileName= resFileName, dirOutput = dirOutput)

# resFileName="simulation_unimodal/parametric/res/2PL_stan_beta.rds"
# dirOutput="simulation_unimodal/parametric/elaboratedRes/"
# args <- list(resFileName= resFileName, dirOutput = dirOutput)
#######
print(str(args))

## results directory
resFileName <- args$resFileName
## output directory
dirOutput <- args$dirOutput
## parametric flag
parametricFlag <- grepl("parametric", resFileName)

modelType <- strsplit(resFileName, "\\/|.rds")[[1]][2]
outFile <- paste0(dirOutput, "efficiency/Efficiency_res_", modelType, ".txt")
################################################################
resObj <- readRDS(resFileName)

tmp <-  strsplit(basename(resFileName), ".rds")[[1]]
label <- strsplit(tmp, "2PL_")[[1]][2]
label <- gsub("_bnp", "", label)

rescale <- TRUE

if(grepl("stan", label)){  
  modelRes <- posteriorRescalingBeta(samples = resObj$samples[, -grep("^eta", colnames(resObj$samples))],
                                   samples2 = resObj$samples[, grep("^eta", colnames(resObj$samples))],
                                   thinEta = 1, 
                                   rescale = rescale)

} else if(grepl("beta", label)) {

	if(grepl("constrained_item", label)) { rescale <- FALSE}

    modelRes <- posteriorRescalingBeta(samples = resObj$samples$samples,
                             samples2 = resObj$samples$samples2,
                             thinEta = resObj$MCMCcontrol$thin2, 
                             rescale = rescale)

  	label <- gsub("_beta", '', label)
  	label <- paste0("beta_", label)

} else if(grepl("gamma", label)){
  
  vars <- colnames(resObj$samples$samples)[grep("lambda", colnames(resObj$samples$samples))]
  
  if(grepl("constrained_item", label)) { rescale <- FALSE }

  label <- gsub("_gamma", '', label)
  label <- paste0("gamma_", label)
  
  modelRes <- posteriorRescalingGamma(samples = resObj$samples$samples,
                             samples2 = resObj$samples$samples2,
                             thinEta = resObj$MCMCcontrol$thin2, 
                             rescale = rescale)
  
  modelRes$betaSamp <- -modelRes$gammaSamp/modelRes$lambdaSamp
  modelRes$betaSamp <-  modelRes$betaSamp - apply(modelRes$betaSamp, 1, mean) 
  colnames(modelRes$betaSamp) <- gsub("lambda", "beta", vars)  
}

saveRDS(modelRes, file = paste0(dirOutput, "inference/", label, "_", modelType,  "_posteriorSamples.rds"))
#####################################################
## Efficency results 
# init values

ess_coda <- NA
ess_mcmcse <- NA
ess_multi <- NA

xx <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
            modelRes[grepl("beta", names(modelRes))][[1]])

## compute ess for item parameters
ess_coda <- min(coda::effectiveSize(xx))
# ess_mcmcse <- min(mcmcse::ess(xx))
ess_multi <- mcmcse::multiESS(xx) ## NaN negative eigenvalue

# ess_coda_all <- coda::effectiveSize(posteriorSamples)
# ess_mcmcse_all <- mcmcse::ess(posteriorSamples)
# ess_multi_all <- mcmcse::multiESS(posteriorSamples) ## NaN negative eigenvalue

path <- paste0(dirOutput, "efficiency/", label, "_", modelType)
compilationTime <- resObj$compilationTime[3]
runningTime <- resObj$runningTime[3]
samplingTime <- 0


## if parametric save also sampling time
if(parametricFlag){ 
  if(grepl("stan", label)){ 
    samplingTime <- resObj$samplingTime
  } else {
    percBurnin <- resObj$MCMCcontrol$nburnin/resObj$MCMCcontrol$niter
    samplingTime <- runningTime * (1 - percBurnin)

  }

}

row <- cbind(path, label, ess_coda, ess_mcmcse, ess_multi, compilationTime, runningTime, samplingTime)

if(!file.exists(outFile)){
	cat(colnames(row), "\n", file = outFile)
}
# append row
cat(row, "\n", file = outFile, append = TRUE)
# ############################################################################