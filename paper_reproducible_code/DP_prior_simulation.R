# library(nimble)
expectedNumberOfClusters <- function(N, alpha){
	sum(alpha/(alpha + seq(1:N) - 1))
}

varianceNumberOfClusters <- function(N, alpha){
	num <- alpha*(seq(1:N) -1)
	den <- (alpha - 1 + seq(1:N))^2
	sum(num/den)
}

#########################
# N = 14525
N = 7377
# N = 2000
R = 10^4

# alphaVals <- rgamma(R, 2, 4)
# alphaVals <- rgamma(R, 1, 3)
alphaVals <- rgamma(R, 1, 1)
plot(density(alphaVals), xlim = c(0,10))
mean(alphaVals)
var(alphaVals)

kVals <- sapply(alphaVals, function(x) expectedNumberOfClusters(N, x))
kValsVar <- sapply(alphaVals, function(x) varianceNumberOfClusters(N, x))

mean(kVals)
## law of total variance! Mean(var(x | alpha)) + var(mean(x|alpha))
mean(kValsVar) + var(kVals)
# mean(kValsVar)

# plot(alphaVals, kVals, col = 2)
# points(alphaVals, kVals - sqrt(kValsVar))
# points(alphaVals, kVals + sqrt(kValsVar))

# plot(density(kVals), xlim = c(0, 10))
# lines(density(kVals[10001:20001]))
# lines(density(kVals[20001:30001]))
# lines(density(kVals[40001:50001]))

expectedNumberOfClusters(N, 1)
varianceNumberOfClusters(N, 1)
log(N)

expectedNumberOfClusters(N, 0.3)
varianceNumberOfClusters(N, 0.3)
log(N)*0.4

expectedNumberOfClusters(10,  0.3)
expectedNumberOfClusters(N,  0.3)
varianceNumberOfClusters(N,  0.3)


##########
alphaVals <- rgamma(R, 2, 6)
plot(density(alphaVals), xlim = c(0,10))
kValsEW <- sapply(alphaVals, function(x) expectedNumberOfClusters(N, x))
kValsVarEW <- sapply(alphaVals, function(x) varianceNumberOfClusters(N, x))

plot(density(kVals), xlim = c(0, 15), main = "Distribution of the expected number of clusters induced by alpha")

lines(density(kValsEW), col = 2)

mean(kValsEW)
mean(kValsVarEW)
# mean(kValsVar)
######################
# Escobar & West -- alpha ~ gamma(2, 4)
# E(alpha) = 0.5; var(alpha) = 0.12

# alpha ~ gamma(1, 3)
# E(alpha) = 0.3; var(alpha) = 0.11
######################
