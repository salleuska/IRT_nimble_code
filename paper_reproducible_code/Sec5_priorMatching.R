################################################
## Function that simulate from the nimble model and return 
## of the marginal prior probability of a correct response
################################################
simulate_samples <- function(model) {
  model$simulate()
  return(model$pi)
}
################################################
library(nimble)

################################################
## Parametric models
################################################
## Hyparameter values are already set in the model scripts for item parameters

parametricModels <- list.files("models/parametric/", full.name = TRUE)

for(i in 1:length(parametricModels)){
	modelName <- parametricModels[[i]]

	## directory for output
	outDir <- "prior_simulations/simulated_samples/"

	modelType       <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[1]
	modelParam      <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[2]
	modelConstraint <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[3]

	outFileName <- paste0(outDir, modelType, "/priorSamples_", modelParam, "_", modelConstraint)


	## if does not exists create directory for output
	dir.create(file.path(outDir, modelType), showWarnings = FALSE)

	## tryCatch to ignore the fact that data is missing and reuse the model code in "models" folder
	tryCatch(source(modelName), 
	error = function(e) message("Ignore error  ", as.character(e)))

	constants <- list(I= 10, N = 100)
	
	model <- nimbleModel(code2PL, constants = constants)
	c_model <- compileNimble(model)
    
	set.seed(20201021 + i)
	t <- system.time(samps <- replicate(100, simulate_samples(c_model)))

	saveRDS(samps, file = paste0(outFileName, ".rds"))	
}

################################################
## BNP models
################################################
## Hyparameter values are already set in the model scripts for item parameters
## In this code we control for hyperpriors for the Dirichlet Process mixture model 

bnpModels <- list.files("models/bnp/", full.name = TRUE)

for(i in 1:length(bnpModels)){
	modelName <- bnpModels[[i]]

	## directory for output
	outDir <- "prior_simulations/simulated_samples/"

	modelType       <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[1]
	modelParam      <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[2]
	modelConstraint <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[3]

	outFileName <- paste0(outDir, modelType, "/priorSamples_", modelParam, "_", modelConstraint)


	## if does not exists create directory for output
	dir.create(file.path(outDir, modelType), showWarnings = FALSE)

	## tryCatch to ignore the fact that data is missing and reuse the model code in "models" folder

	constants <- list(I= 10, N = 100, M =100)

	tryCatch(source(modelName), 
	error = function(e) message("Ignore error  ", as.character(e)))

	##-----------------------------------------##
	## Dirichlet process mixture initialization
	##-----------------------------------------##
	inits <- list(nu1 = 2.01, nu2 = 1.01, ## s2 ~ InvGamma(nu1, nu2)
				  s2_mu = 2,              ## mu ~ N(0, s2_mu)
	              a = 2, b = 4)           ## Escobar & West prior for DP conc   

	model <- nimbleModel(code2PL, 
						 constants = constants, 
						 inits = inits,
						 calculate = FALSE)

	c_model <- compileNimble(model)
    
	set.seed(20201021 + i)
	t <- system.time(samps <- replicate(100, simulate_samples(c_model)))

	saveRDS(samps, file = paste0(outFileName, "_a_",c_model$a, "_b_",c_model$b,".rds"))	

	##-----------------------------------------##
	## other prior for DP conc
	##-----------------------------------------##

	c_model$a <- 1
	c_model$b <- 3

	set.seed(20201022 + i)
	t <- system.time(samps <- replicate(100, simulate_samples(c_model)))

	saveRDS(samps, file = paste0(outFileName, "_a_",c_model$a, "_b_",c_model$b,".rds"))	

	##-----------------------------------------##
	## other prior for DP conc
	##-----------------------------------------##
	c_model$a <- 1
	c_model$b <- 1

	set.seed(20201022 + i)
	t <- system.time(samps <- replicate(100, simulate_samples(c_model)))

	saveRDS(samps, file = paste0(outFileName, "_a_",c_model$a, "_b_",c_model$b,".rds"))	
}



################################################
## BNP models - fixed alpha
################################################
## Hyparameter values are already set in the model scripts for item parameters
## In this code we control different fixed values of the DP conc parameters

bnpModels <- list.files("models/bnp_fixedAlpha/", full.name = TRUE)
alphaVec <- c(0.01, 0.05, 0.5, 1, 1.5, 2)

for(i in 1:length(bnpModels)){
	modelName <- bnpModels[[i]]

	## directory for output
	outDir <- "prior_simulations/simulated_samples/"

	modelType       <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[1]
	modelParam      <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[2]
	modelConstraint <- unlist(strsplit(basename(modelName), "[\\_\\.]"))[3]

	outFileName <- paste0(outDir, modelType, "/priorSamples_", modelParam, "_", modelConstraint)


	## if does not exists create directory for output
	dir.create(file.path(outDir, modelType), showWarnings = FALSE)

	## tryCatch to ignore the fact that data is missing and reuse the model code in "models" folder

	constants <- list(I= 10, N = 100, M =100)

	tryCatch(source(modelName), 
	error = function(e) message("Ignore error  ", as.character(e)))

	##-----------------------------------------##
	## Dirichlet process mixture initialization - fixed alpha
	##-----------------------------------------##
	for(j in 1:length(alphaVec)){
		inits <- list(nu1 = 2.01, nu2 = 1.01, ## s2 ~ InvGamma(nu1, nu2)
					  s2_mu = 2)              ## mu ~ N(0, s2_mu)
		              
		model <- nimbleModel(code2PL, 
							 constants = constants, 
							 inits = inits,
							 calculate = FALSE)
		
		model[["alpha"]] <- alphaVec[j]
		c_model <- compileNimble(model)

		set.seed(20201027 + i + j)
		t <- system.time(samps <- replicate(100, simulate_samples(c_model)))

		saveRDS(samps, file = paste0(outFileName, "_alhpaFixed_", c_model$alpha, ".rds"))			
	}

}
