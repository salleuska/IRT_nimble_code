##----------------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##----------------------------------------------#
## this code run stan model - default IRT 2pl found at
## https://mc-stan.org/users/documentation/case-studies/tutorial_twopl.html#estimate2pl
##--------------------------------------------------------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --model=
## --dirResults=
## --data=
## --niter=
## --nburnin=

data="data/simulation_unimodal.rds" 
niter=2000 
nburnin=2000 

args <- list(data = data, niter = niter, nburnin = nburnin)
##----------------------------------##
## load library and functions
library(rstan)
library(reshape2)

##----------------------------------##
cat("Warning: using burning as n. warmup iterations. \n
	N. iterations will be nburnin + niter")

## Load data
if(grepl("timss", args$data)){
	alldata <- readRDS(args$data)
	data <- list(y = alldata$y)
} else {
	data <- list(y = readRDS(args$data))
}

##-----------------------------------------##
## Set variables 
##-----------------------------------------##
## results directory
if(is.null(args$dirResults)) dir <- "output/posterior_samples" else dir <- args$dirResults

MCMCcontrol <- list()

MCMCcontrol$nwarmup <- as.numeric(args$nburnin)
MCMCcontrol$niter <- as.numeric(args$niter) + MCMCcontrol$nwarmup

fileStan <- "models/parametric_IRT_constrainedAbilities.stan"

##----------------------------------##
cat("Stan 2PL IRT constrained abilities model\n")
##----------------------------------##

if(grepl("timss", args$data)){

	stan_data <- list(I  = max(alldata$item), 
					  J  = max(alldata$id),
					  N  = nrow(alldata), 
					  ii = alldata$item, 
	                  jj = alldata$id, 
	                  y  = alldata$y)

} else {

	## Reshape data in wide format
	wide    <- as.data.frame(data)
	wide$id <- 1:nrow(wide)  # Attach a person ID number to each row.
	long    <- melt(wide, id.vars = "id", variable.name = "item", value.name = "response")
	# head(long)

	key <- 1:length(unique(long$item))
	names(key) <- unique(long$item)
	long$item.id <- key[long$item]
	# head(long)

	stan_data <- list(I  = max(long$item.id), 
					  J  = max(long$id),
					  N  = nrow(long), 
					  ii = long$item.id, 
	                  jj = long$id,
	                  y  = long$response)
}

## Create arguments lists
stan_model_args <- list() 
sampling_args   <- list()

## modify stan_model_args
stan_model_args$file <- fileStan
## Create stan_model object
compileTime <- system.time(stan_mod <- do.call(rstan::stan_model, stan_model_args))

## modify sampling args
sampling_args$object <- stan_mod ## object of class stanmodel
sampling_args$data   <- stan_data
sampling_args$chains <- 1

##  Note: in rstan::sampling function the `iter` argument comprises also the number of warmup iterations
sampling_args$warmup <- MCMCcontrol$nwarmup
sampling_args$iter   <- MCMCcontrol$niter 
sampling_args$thin   <- 1
sampling_args$seed   <- 1
     
totalTime <- system.time(stan_out <- do.call(rstan::sampling, sampling_args))

monitors <- c("beta", "lambda", "eta")

samplesArray <- rstan::extract(stan_out, 
                            pars = monitors,
                            permuted = FALSE,
                            inc_warmup = FALSE)[, 1, ]

sampleTime <- rstan::get_elapsed_time(stan_out)[2]

results <- list(samples          = samplesArray,
                compilationTime  = compileTime,
                runningTime      = totalTime,
                samplingTime     = sampleTime,
                MCMCcontrol      = MCMCcontrol)

## directory for output
modelType       <- unlist(strsplit(basename(fileStan), "[\\_\\.]"))[1]
dataName        <- unlist(strsplit(basename(args$data), "[.]"))[1]

outDir <- paste0(dir, "/", dataName, "/", modelType, "/")

dir.create(file.path(outDir), recursive = TRUE, showWarnings = FALSE)

saveRDS(results, file  = paste0(outDir, "parametric_IRT_stan.rds"))
