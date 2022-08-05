## ---- simulation1
## Simulated data with unimodal latent abilities
set.seed(1)

# set the number of items I and persons P
nItems    <- 15
nPersons  <- 2000

# set the item and population parameters
lambda0     <- c(1 + runif(nItems,-0.5,0.5))
beta0     <- seq(-3,3,length=nItems)


lambda0 <- exp(log(lambda0) - mean(log(lambda0)))

etaMean    <- 0
etaSigma2 <- (1.25)^2

# generate thetas and the I x P matrix of response probabilities
etaAbility <- rnorm(nPersons, mean=etaMean, sd=sqrt(etaSigma2))

term1     <- outer(etaAbility, lambda0)
term2     <- matrix(rep(lambda0*beta0, nPersons), nrow=nPersons,byrow=TRUE)
prob      <- plogis(term1-term2)  ### 1/(1 + exp(term.2 - term.1))
# generate the 0/1 responses Y as a matrix of Bernoulli draws

Y        <- ifelse(runif(nItems*nPersons) < prob, 1, 0)
save(nItems, nPersons, lambda0, beta0, etaAbility, Y , prob, file = "simulation_unimodal_allValues.RData")
saveRDS(Y , file = "simulation_unimodal.rds")
## ---- end-of-simulation1

## ---- simulation1 Expanded
## Simulated data with unimodal latent abilities
set.seed(243)

# set the number of items I and persons P
nItems    <- c(10, 30)
nPersons  <- c(1000, 5000)

for(i in 1:2){
	for(j in 1:2) {
	
	nItemSim = nItems[i]
	nPersonSim = nPersons[j]

	# set the item and population parameters
	lambda0     <- c(1 + runif(nItemSim,-0.5,0.5))
	beta0     <- seq(-3,3,length=nItemSim)


	lambda0 <- exp(log(lambda0) - mean(log(lambda0)))

	etaMean    <- 0
	etaSigma2 <- (1.25)^2

	# generate thetas and the I x P matrix of response probabilities
	etaAbility <- rnorm(nPersonSim, mean=etaMean, sd=sqrt(etaSigma2))

	term1     <- outer(etaAbility, lambda0)
	term2     <- matrix(rep(lambda0*beta0, nPersonSim), nrow=nPersonSim,byrow=TRUE)
	prob      <- plogis(term1-term2)  ### 1/(1 + exp(term.2 - term.1))
	# generate the 0/1 responses Y as a matrix of Bernoulli draws

	Y        <- ifelse(runif(nItemSim*nPersonSim) < prob, 1, 0)
	
	filename <- paste0("simulation_unimodal_I_", nItemSim, "_N_", nPersonSim)
	save(nItemSim, nPersonSim, lambda0, beta0, etaAbility, Y , prob, file = paste0(filename, "_allData.RData"))
	saveRDS(Y , file = paste0(filename, ".rds"))


	}
}


## ---- end-of-simulation1-expanded
