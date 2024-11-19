##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
library(sn)

dMultiModal <- function(x, weights = c(1,1,1), means = c(-2, 0, 3)){
    require(sn)

    prop <- weights/sum(weights)

    prop[1]*dnorm(x, mean = means[1], sd = sqrt(1)) + 
    prop[2]*dnorm(x, mean = means[2], sd = sqrt(0.5)) + 
    prop[3]*dsn(x, xi = means[3], omega = 1, alpha = -3)
}

pMultiModal <- function(x, weights = c(1,1,1), means = c(-2, 0, 3)){
    require(sn)

    prop <- weights/sum(weights)

    prop[1]*pnorm(x, mean = means[1], sd = sqrt(1), lower.tail=TRUE) + 
    prop[2]*pnorm(x, mean = means[2], sd = sqrt(0.5), lower.tail=TRUE) + 
    prop[3]*psn(x, xi = means[3], omega = 1, alpha = -3)
}



rMultiModal <- function(n, weights = c(1,1,1), means = c(-2, 0, 3)){
    prop <- weights/sum(weights)

    out <- c(rnorm(n*prop[1], mean = means[1], sd = sqrt(1)),
    rnorm(n*prop[2], mean = means[2], sd = sqrt(0.5)),
    rsn(n*prop[3], xi = means[3], omega = 1, alpha = -3))

    out
}
