##----------------------------##
## Multimodal
##----------------------------##
library(sn)

dMultiModal <- function(x, weights = c(1,1,1), means = c(-2, 0, 3)){
	prop <- weights/sum(weights)

	prop[1]*dnorm(x, mean = means[1], sd = sqrt(1)) + 
	prop[2]*dnorm(x, mean = means[2], sd = sqrt(0.5)) + 
	prop[3]*dsn(x, xi = means[3], omega = 1, alpha = -3)
}

## curve(dMultiModal(x, weights = c(4,4,2)), from = -5, to = 5)

rMultiModal <- function(n, weights = c(4,4,2), means = c(-2, 0, 3)){
	prop <- weights/sum(weights)

	out <- c(rnorm(n*prop[1], mean = means[1], sd = sqrt(1)),
	rnorm(n*prop[2], mean = means[2], sd = sqrt(0.5)),
	rsn(n*prop[3], xi = means[3], omega = 1, alpha = -3))

	out
}
## hist(rMultiModal(n = 10000), breaks = 100)

## ---- simulation
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

## parameters of the 3 components mixture

# possible change
weights = c(2,4,4)
means = c(-2, 0, 3)

# generate thetas from a mixture and the I x P matrix of response probabilities
etaAbility  <- rMultiModal(n = nPersons, weights = weights, means = means)
plot(density(etaAbility))

term1       <- outer(etaAbility, lambda0)
term2       <- matrix(rep(lambda0*beta0, nPersons), nrow=nPersons,byrow=TRUE)
prob        <- plogis(term1-term2)  ### 1/(1 + exp(term.2 - term.1))
# generate the 0/1 responses Y as a matrix of Bernoulli draws

Y        <- ifelse(runif(nItems*nPersons) < prob, 1, 0)
save(nItems, nPersons, lambda0, beta0, etaAbility, Y , prob, file = "simulation_multimodal2_allValues.RData")
saveRDS(Y , file = "simulation_multimodal2.rds")


## ---- simulation Expanded
## Simulated data with multimodal2 latent abilities
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
	etaAbility  <- rMultiModal(n = nPersonSim, weights = weights, means = means)

	term1     <- outer(etaAbility, lambda0)
	term2     <- matrix(rep(lambda0*beta0, nPersonSim), nrow=nPersonSim,byrow=TRUE)
	prob      <- plogis(term1-term2)  ### 1/(1 + exp(term.2 - term.1))
	# generate the 0/1 responses Y as a matrix of Bernoulli draws

	Y        <- ifelse(runif(nItemSim*nPersonSim) < prob, 1, 0)
	
	filename <- paste0("simulation_multimodal2_I_", nItemSim, "_N_", nPersonSim)
	save(nItemSim, nPersonSim, lambda0, beta0, etaAbility, Y , prob, file = paste0(filename, "_allData.RData"))
	saveRDS(Y , file = paste0(filename, ".rds"))


	}
}


## ---- end-of-simulation1-expanded

