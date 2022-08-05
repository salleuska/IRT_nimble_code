##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
# This files compute some quantities for plots and tables
dir.create("figures/dataForFigures", recursive = TRUE, showWarnings = FALSE)
##-----------------------------------------#
## Data
# rm(list = ls()# 
library(bayestestR) ## For HDI intervals
# ## Compute simulation MSE
source("R_functions/ggplot_settings.R")


##--------------------------------##
## Unimodal simulation
##--------------------------------##
## load true values
dataName <- "simulation_unimodal"
load(paste0("data/",dataName,"_allValues.RData"))
## Set a grid for density computation
grid <- seq(-8, 8, len = 400) 

trueValues <- list(beta   = beta0,
		   		   lambda = lambda0,
				     eta    = etaAbility)

bestModel <- modelData$model[modelData$data == dataName]

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


biasParametric <-  sapply(1:3, function(i) paraEstimates[[i]] - trueValues[[i]])
biasBnp        <-  sapply(1:3, function(i) bnpEstimates[[i]] - trueValues[[i]])

metricsUnimodal <- data.frame(parameters =c("Difficulties", "Discrimination", "Abilities"),  
							   MAE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(abs(x)))), 
							   MSE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(x^2))),
							   MAE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(abs(x)))), 
							   MSE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(x^2))))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))

##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)
 
## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
indices <- seq(10, 45000, by = 10)
bnpG0 <- bnpG0[indices]

densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3])))*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
truePerc <- pnorm(etaAbility , 0, sd = 1.25)

indexSample <- order(truePerc)[round(seq(1, 2000, length = 50))]
## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]

## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}
##########################

unimodalRes <- list(truValues = trueValues,
					paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure, 
					truePerc = truePerc[indexSample],
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(unimodalRes, file = "figures/dataForFigures/unimodal.rds")

##--------------------------------##
## Bimodal simulation
##--------------------------------##
## load true values
dataName <- "simulation_bimodal"
load(paste0("data/",dataName,"_allValues.RData"))

## Set a grid for density computation
grid <- seq(-8, 8, len = 400) 

trueValues <- list(beta   = beta0,
		   		   lambda = lambda0,
				   eta    = etaAbility)

bestModel <- modelData$model[modelData$data == dataName]

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


biasParametric <-  sapply(1:3, function(i) paraEstimates[[i]] - trueValues[[i]])
biasBnp        <-  sapply(1:3, function(i) bnpEstimates[[i]] - trueValues[[i]])

## Save simulation metrics computed using true values
metricsBimodal <- data.frame(parameters =c("Difficulties", "Discrimination", "Abilities"),  
							   MAE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(abs(x)))), 
							   MSE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(x^2))),
							   MAE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(abs(x)))), 
							   MSE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(x^2))))

save(metricsUnimodal, metricsBimodal, file = "figures/dataForFigures/table3_metricsSimulations.RData")
#--------------------------------------------------------------------------------#
# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2,function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#

niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid + paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
indices <- seq(10, 45000, by = 10)
bnpG0 <- bnpG0[indices]
densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute percentiles simulation - qqplot
##------------------------------------------------------------#
## Take a grid
##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
truePerc <- 0.5*pnorm(etaAbility , -2, sd = 1.25) + 0.5*pnorm(etaAbility , 2, sd = 1.25)

indexSample <- order(truePerc)[round(seq(1, 2000, length = 50))]
## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]


## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}


# plot(truePerc[indexSample])
# points(1:50, apply(paraPerc, 2, mean), col = 4)
# points(1:50, apply(bnpPerc, 2, mean), col = 3)
##-----------------------------------------------##

bimodalRes <- list(truValues = trueValues,
					paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					truePerc = truePerc[indexSample],
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(bimodalRes, file = "figures/dataForFigures/bimodal.rds")
##------------------------------------------------------------#
##--------------------------------##
## multimodal simulation
##--------------------------------##
## load true values
dataName <- "simulation_multimodal"
#dataName <- "simulation_multimodal"
load(paste0("data/",dataName,"_allValues.RData"))
## Set a grid for density computation
grid <- seq(-8, 8, len = 400) 

trueValues <- list(beta   = beta0,
		   		   lambda = lambda0,
				     eta    = etaAbility)

bestModel <- modelData$model[modelData$data == dataName]


bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


biasParametric <-  sapply(1:3, function(i) paraEstimates[[i]] - trueValues[[i]])
biasBnp        <-  sapply(1:3, function(i) bnpEstimates[[i]] - trueValues[[i]])

metricsMultimodal <- data.frame(parameters =c("Difficulties", "Discrimination", "Abilities"),  
							   MAE_multimodal_para = unlist(lapply(biasParametric, function(x) mean(abs(x)))), 
							   MSE_multimodal_para = unlist(lapply(biasParametric, function(x) mean(x^2))),
							   MAE_multimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(abs(x)))), 
							   MSE_multimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(x^2))))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)
# 
## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
indices <- seq(10, 45000, by = 10)
bnpG0 <- bnpG0[indices]

densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3])))*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
truePerc <- pnorm(etaAbility , 0, sd = 1.25)

indexSample <- order(truePerc)[round(seq(1, 2000, length = 50))]
## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]

## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}
##########################

multimodalRes <- list(truValues = trueValues,
					paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure, 
					truePerc = truePerc[indexSample],
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(multimodalRes, file = "figures/dataForFigures/multimodal.rds")
##------------------------------------------------------------#
## Data
##------------------------------------------------------------#
##--------------------------------##
## Health data
##--------------------------------##
rm(list = ls());gc()
library(bayestestR) ## For HDI intervals

## Compute simulation MSE
source("R_functions/ggplot_settings.R")

## Set a grid for density computation
grid <- seq(-15, 30, len = 800) 

dataName <- "data_health"

bestModel <- modelData$model[modelData$data == dataName]
# bestModel <- "SI_unconstrained"

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#

niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  # rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
niter2 <- length(bnpG0)
indices <- seq(10, niter2, by = 10)

bnpG0 <- bnpG0[indices]
densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
	# rescGrid <- bnpModel$scaleShiftEta[i]*grid -  bnpModel$locationShiftEta[i]
	rescGrid <- bnpModel$scaleShiftEta[i]*grid +  bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
indexSample <- order(bnpEstimates$eta)[round(seq(1, length(bnpEstimates$eta), length = 50))]

## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]

## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

#####

healthRes   <- list(paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)


saveRDS(healthRes, file = "figures/dataForFigures/health.rds")

##--------------------------------##
## Timss data
##--------------------------------##
rm(list = ls());gc()
library(bayestestR) ## For HDI intervals

## Compute simulation MSE
source("R_functions/ggplot_settings.R")

## Set a grid for density computation
grid <- seq(-8, 8, len = 800) 

dataName <- "data_timss"
bestModel <- modelData$model[modelData$data == dataName]

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
niter2 <- length(bnpG0)
indices <- seq(10, niter2, by = 10)

bnpG0 <- bnpG0[indices]

densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid -  bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}
##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#


indexSample <- order(bnpEstimates$eta)[round(seq(1, length(bnpEstimates$eta), length = 50))]
## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]


## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] - paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] - bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

##------------------------------------------------------------#

timssRes   <- list(paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)


saveRDS(timssRes, file = "figures/dataForFigures/timss.rds")
##------------------------------------------------------------#
## TIMSS 3PLx
##------------------------------------------------------------#
rm(list = ls());gc()
library(bayestestR) ## For HDI intervals

## Compute simulation MSE
source("R_functions/ggplot_settings.R")

## Set a grid for density computation
grid <- seq(-10, 8, len = 800) 

dataName <- "data_timss"
bestModel <- "IRT_unconstrained"

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp3PL/bnp3PL_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric3PL/parametric3PL_", bestModel, ".rds"))

paraModel$deltaSamp <- paraModel$otherParSamp[, grep("delta", colnames(paraModel$otherParSamp))]
bnpModel$deltaSamp  <- bnpModel$otherParSamp[, grep("delta", colnames(paraModel$otherParSamp))]


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
                  lambda = apply(paraModel$lambdaSamp, 2, mean), 
                  eta    = apply(paraModel$etaSamp, 2, mean), 
                  delta  = apply(paraModel$deltaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
                      lambda = apply(bnpModel$lambdaSamp, 2, mean), 
                      eta    = apply(bnpModel$etaSamp, 2, mean), 
                      delta  = apply(bnpModel$deltaSamp, 2, mean))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
                lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
                eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low), 
                delta    = apply(paraModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))

paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
                   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
                    eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high),
                    delta   = apply(paraModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


bnpLow <-   list(beta    = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
                lambda   = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
                eta      = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low), 
                delta    = apply(bnpModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
                   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
                   eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
                   delta  = apply(bnpModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))

##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid + paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp3PL_", bestModel, ".rds"))
niter2 <- length(bnpG0)

indices <- seq(10, niter2, by = 10)

bnpG0 <- bnpG0[indices]
densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
        rescGrid <- bnpModel$scaleShiftEta[i]*grid +  bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}
##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
indexSample <- order(bnpEstimates$eta)[round(seq(1, length(bnpEstimates$eta), length = 50))]
## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]


## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
         rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
         paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
         rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
         bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

##------------------------------------------------------------#
timssRes   <- list(paraEstimates = paraEstimates, 
                                        paraLow = paraLow, 
                                        paraUpper = paraUpper,
                                        bnpEstimates = bnpEstimates, 
                                        bnpLow = bnpLow, 
                                        bnpUpper = bnpUpper,
                                        grid = grid, 
                                        densitySamplesPara = densitySamplesPara, 
                                        densityDPMeasure = densityDPMeasure,
                                        paraPerc = paraPerc, 
                                        bnpPerc = bnpPerc)

saveRDS(timssRes, file = "figures/dataForFigures/timss3PL.rds")