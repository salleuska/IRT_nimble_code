##################################################
## These script containts functions to compute the expectation and variance
## of the expected number of clusters via Monte Carlo when using a prior 
## distribution for the DP concentration parameter alpha
##################################################

expectedNumberOfClusters <- function(N, alpha){
	sum(alpha/(alpha + seq(1:N) - 1))
}

varianceNumberOfClusters <- function(N, alpha){
	num <- alpha*(seq(1:N) -1)
	den <- (alpha - 1 + seq(1:N))^2
	sum(num/den)
}

#########################
## Example 
#########################
## numer of observations in the data
N <- 2000
## parameters for the gamma distribution
a <- 1
b <- 1

## number of Monte carlo replicates
R <- 10^5

## generate alpha ~ Ga(a, b)
set.seed(1)
alphaVals <- rgamma(R, a, b)

kVals <- sapply(alphaVals, function(x) expectedNumberOfClusters(N, x))
kValsVar <- sapply(alphaVals, function(x) varianceNumberOfClusters(N, x))

## Estimate of the a priori expected number of clusters
mean(kVals)

## Estimate of the a priori variance for thenumber of clusters
## uses law of total variance! Mean(var(x | alpha)) + var(mean(x|alpha))
mean(kValsVar) + var(kVals)

