##----------------------------##
## mixture model
##----------------------------##
## ---- simulation2
# Set the random-number generator seed,
# to make the results reproducible
set.seed(10)

# set the number of items I and persons P
nItems    <- 15
nPersons  <- 2000

# set the item and population parameters
lambda0     <- c(1 + runif(nItems,-0.5,0.5))
beta0       <- seq(-3,3,length=nItems)


lambda0 <- exp(log(lambda0) - mean(log(lambda0)))

etaMean1  <- -2
etaMean2 <- 2
etaSigma2 <- (1.25)^2


# generate thetas from a mixture and the I x P matrix of response probabilities
etaAbility  <- rnorm(nPersons/2, mean=etaMean1, sd=sqrt(etaSigma2))
etaAbility  <- c(etaAbility, rnorm(nPersons/2, mean=etaMean2, sd=sqrt(etaSigma2)))
term1       <- outer(etaAbility, lambda0)
term2       <- matrix(rep(lambda0*beta0, nPersons), nrow=nPersons,byrow=TRUE)
prob        <- plogis(term1-term2)  ### 1/(1 + exp(term.2 - term.1))
# generate the 0/1 responses Y as a matrix of Bernoulli draws

Y        <- ifelse(runif(nItems*nPersons) < prob, 1, 0)
save(nItems, nPersons, lambda0, beta0, etaAbility, Y , prob, file = "simulation_bimodal_allValues.RData")
saveRDS(Y , file = "simulation_bimodal.rds")
## ---- end-of-simulation2

## ---- simulation2 Expanded
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

	# generate thetas from a mixture and the I x P matrix of response probabilities
	etaAbility  <- rnorm(nPersonSim/2, mean=etaMean1, sd=sqrt(etaSigma2))
	etaAbility  <- c(etaAbility, rnorm(nPersonSim/2, mean=etaMean2, sd=sqrt(etaSigma2)))
	term1       <- outer(etaAbility, lambda0)
	term2       <- matrix(rep(lambda0*beta0, nPersonSim), nrow=nPersonSim,byrow=TRUE)
	prob        <- plogis(term1-term2)  ### 1/(1 + exp(term.2 - term.1))

	# generate the 0/1 responses Y as a matrix of Bernoulli draws

	Y        <- ifelse(runif(nItemSim*nPersonSim) < prob, 1, 0)
	
	filename <- paste0("simulation_bimodal_I_", nItemSim, "_N_", nPersonSim)
	save(nItemSim, nPersonSim, lambda0, beta0, etaAbility, Y , prob, file = paste0(filename, "_allData.RData"))
	saveRDS(Y , file = paste0(filename, ".rds"))


	}
}

