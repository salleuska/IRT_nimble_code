##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
## This scripts reconstruct a BNP model & relative MCMC
## populate the MCMC with posterior samples, and 
## use getSamplesDPMeasure to sample from the real DP
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --dataName
## --modelName

## dataName <- "simulation_bimodal"
## modelName <- "bnp_IRT_unconstrained"

## Set directories
dirResults <- "output/posterior_samples/"
dirOutput  <- "output/posterior_samples_elaborated/"

## dataName
dataName <- args$dataName
## model
modelName <- args$modelName
## Path to results 
path <- paste0(dirResults, dataName, "/bnp/")
##-----------------------------------------#
## Load libraries
library(nimble)
##-----------------------------------------#
## Read data
if(grepl("timss", dataName)){
	alldata <- readRDS(paste0("data/",dataName,".rds"))
	data <- list(y = alldata$y)
} else {
	data <- list(y = readRDS(paste0("data/",dataName,".rds")))
}

## Read original results 
originalRes     <- readRDS(paste0(path, modelName, ".rds"))
originalSamples <- originalRes$samples[[1]]

## thinning for eta
thinEta <- originalRes$MCMCcontrol$thin2
nSamp <- originalRes$MCMCcontrol$niter - originalRes$MCMCcontrol$nburnin

indicesEta <- seq(from = thinEta,  
				  to = nSamp, 
				  by = thinEta)

## use iterations corresponding to r.e. thinning
samples <- originalSamples[indicesEta, ]

## reconstruct MCMC samples matrix for the bnp model
samplesForDP <- samples[, grep(c("alpha|zi|muTilde|s2Tilde"), colnames(samples))]

##-----------------------------------------#
## Read model
##-----------------------------------------#

if(grepl("timss", dataName)){
	source(paste0("models/bnp_timss/", modelName, ".R"))
	constants$M <- 30
} else {
	source(paste0("models/bnp/", modelName, ".R"))
}

## create and compile model
model  <- nimbleModel(code2PL, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)
## simulate
cmodel$simulate()

conf <- configureMCMC(model)
conf$monitors <- c("alpha", "zi", "muTilde", "s2Tilde")	

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)

## put samples back into compiled mcmcm object
nimble:::matrix2mv(samplesForDP, cmcmc$mvSamples)

##-----------------------------------------#
## Samples from G0 base measure for abilities 
##-----------------------------------------#
t <- system.time(outputG <- getSamplesDPmeasure(cmcmc))

outDir <- paste0(dirOutput, dataName)
dir.create(file.path(outDir), recursive = TRUE, showWarnings = FALSE)

saveRDS(outputG, file = paste0(outDir, "/DPG0_", modelName, ".rds"))
##################################
# Rscript simulateFromDPmeasure.R \
# --dataName="simulation_unimodal" \
# --modelName="2PL_unconstrained_bnp_beta"

# Rscript simulateFromDPmeasure.R \
# --dataName="simulation_bimodal" \
# --modelName="2PL_unconstrained_bnp_beta"

# Rscript simulateFromDPmeasure.R \
# --dataName="simulation_bimodal_skewed" \
# --modelName="2PL_unconstrained_bnp_beta"

# Rscript simulateFromDPmeasure.R \
# --dataName="data_health" \
# --modelName="2PL_unconstrained_bnp_beta"

# Rscript simulateFromDPmeasure.R \
# --dataName="data_health" \
# --modelName="2PL_unconstrained_bnp_gamma"

# Rscript simulateFromDPmeasure.R \
# --dataName="data_health" \
# --modelName="2PL_constrained_bnp_item_beta"

# Rscript simulateFromDPmeasure.R \
# --dataName="data_timss" \
# --modelName="2PL_unconstrained_bnp_beta"

# Rscript simulateFromDPmeasure.R \
# --dataName="data_timss" \
# --modelName="2PL_unconstrained_bnp_gamma"
####
