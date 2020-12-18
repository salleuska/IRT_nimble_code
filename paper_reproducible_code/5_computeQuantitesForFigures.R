##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
## This files compute some quantities for plots and tables
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
indices <- seq(10, 45000, by = 10)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[indices ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[indices ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
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
## Compute percentiles simulation - qqplot
##------------------------------------------------------------#
## Take a grid
gridQuantiles <- qnorm(seq(0, 1, length = 50), 0, 1.25)
paraPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*gridQuantiles - paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}
bnpPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*gridQuantiles - bnpModel$locationShiftEta[i]
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
					gridQuantiles = gridQuantiles,
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
indices <- seq(10, 45000, by = 10)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[indices ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[indices ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid + paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
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

gridQuantiles <- c(-Inf, seq(-5, 5, length = 48), Inf)

paraPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*gridQuantiles + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}
bnpPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*gridQuantiles + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

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
					gridQuantiles = gridQuantiles,
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(bimodalRes, file = "figures/dataForFigures/bimodal.rds")
##------------------------------------------------------------#
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
indices <- seq(10, 45000, by = 10)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[indices ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[indices ,"s2.eta"]

for(i in 1:niter){
  # rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
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
## Compute percentiles simulation - qqplot
##------------------------------------------------------------#
## Take a grid
gridQuantiles <- c(-Inf, seq(-10, 20, length = 48), Inf)
# gridQuantiles <- seq(-10, 20, length = 200)

paraPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*gridQuantiles + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}
bnpPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*gridQuantiles + bnpModel$locationShiftEta[i]
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
					gridQuantiles = gridQuantiles,
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
indices <- seq(10, 45000, by = 10)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[indices ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[indices ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
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
## Compute percentiles simulation - qqplot
##------------------------------------------------------------#
## Take a grid

gridQuantiles <- c(-Inf, seq(-3.5, 3.5, length = 48), Inf)

paraPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*gridQuantiles + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}
bnpPerc  <- matrix(0, ncol = length(gridQuantiles), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*gridQuantiles + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}
#####
##------------------------------------------------------------#
## Compute percentiles simulation - qqplot
##------------------------------------------------------------#
## Take a grid


timssRes   <- list(paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					gridQuantiles = gridQuantiles,
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)


saveRDS(timssRes, file = "figures/dataForFigures/timss.rds")
