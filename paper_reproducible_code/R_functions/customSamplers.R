##------------------------------##
## Centered sampler
##------------------------------##
## Nested sampler implementation

## Define a new base class for individual centered samplers.
## Since they need information about the mean covariate for centering,
## the base class sampler_BASE will not be sufficient.
sampler_centered_single_BASE <- nimbleFunctionVirtual(
    methods = list(
        reset = function() { },
        set_mean = function(m = double()) {}
    )
)

## This is what does the actual sampling operation.
sampler_centered_single <- nimbleFunction(
  name = "centered_single",
  contains = sampler_centered_single_BASE,
  setup = function(model, mvSaved, target, control) {
    centering_mean <- 0
    scale         <- if(!is.null(control$scale))         control$scale         else 1
    ##
    adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
    adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200

 	
 	  ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    if(length(targetAsScalar) != 2)   stop('must use centered sampler on exactly two nodes')

    target1 <- targetAsScalar[1]
    target2 <- targetAsScalar[2]
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)

    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
        saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0

	## not sure if necessary
	if(any(model$isDiscrete(targetAsScalar))) stop('cannot use centered sampler on discrete-valued target')
  },
  run = function() {

     ## This code assumes lambda in original scale
     #    currentValue1 <- model[[target1]]
     #    currentValue2 <- model[[target2]]

    	# propLogScale <- rnorm(1, mean = 0, sd = scale)
     #    propFactor <- exp(propLogScale)
     #    ## lambda proposal
     #    propValue1 <- currentValue1 * propFactor
     #    ## gamma value (centered)
     #    propValue2 <- currentValue2 + centering_mean*(currentValue1 - propValue1)

     #    model[[target1]] <<- propValue1
     #    model[[target2]] <<- propValue2
     #    logMHR <- calculateDiff(model, target)
 
     ## This code assumes lambda in log scale
        currentValue1 <- model[[target1]]
        currentValue2 <- model[[target2]]

    	propLogScale <- rnorm(1, mean = 0, sd = scale)
        ## lambda proposal
        propValue1 <- currentValue1 + propLogScale
        ## gamma value (centered)
        propValue2 <- currentValue2 + centering_mean*(exp(currentValue1) - exp(propValue1))

        model[[target1]] <<- propValue1
        model[[target2]] <<- propValue2
        logMHR <- calculateDiff(model, target)
        if(logMHR == -Inf) {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            ## Drawing a random number is needed during first testing
            ## of this step in order to keep the random numbers identical
            ## to old behavior to see if tests that depend on particular
            ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
            ## we call decide() to ensure same behavior.
            ## jump <- decide(logMHR)
            ## When new behavior is acceptable, we can remove the above line
            ## and uncomment the following:
            jump <- FALSE
        } else {
            ## if lambda in original scale
            # logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf) + propLogScale
           	## if lambda in log scale
            logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf)
            jump <- decide(logMHR)
            if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }  
        if(adaptive)     adaptiveProcedure(jump)

  }, 
  methods = list(
    set_mean = function(m = double()) {
     centering_mean <<- m
    }, 
    adaptiveProcedure = function(jump = logical()) {
    timesRan <<- timesRan + 1
    if(jump)     timesAccepted <<- timesAccepted + 1
    if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
            setSize(scaleHistory, timesAdapted)         ## scaleHistory
            scaleHistory[timesAdapted] <<- scale        ## scaleHistory
            setSize(acceptanceHistory, timesAdapted)         ## scaleHistory
            acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^0.8)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        }

        timesRan <<- 0
        timesAccepted <<- 0
    },
    getScaleHistory = function() {  ## scaleHistory
        returnType(double(1))
        if(saveMCMChistory) {
            return(scaleHistory)
        } else {
            print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
            return(numeric(1, 0))
        }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
        returnType(double(1))
        if(saveMCMChistory) {
            return(acceptanceHistory)
        } else {
            print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
            return(numeric(1, 0))
        }
    },          
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
        scale <<- scaleOriginal
        timesRan      <<- 0
        timesAccepted <<- 0
        timesAdapted  <<- 0
        if(saveMCMChistory) {
            scaleHistory  <<- c(0, 0)    ## scaleHistory
            acceptanceHistory  <<- c(0, 0)
        }
        gamma1 <<- 0
    }) 
)


## This is what will be assigned as a sampler for all the pairs, e.g.
## mcmcConf$addSampler(type = 'centered_samplers', 
##            target = c('lambda', 'gamma'),
##          control = list(nodesToCenter = 'eta') )
## 
sampler_centered <- nimbleFunction(
 name = "centered",
contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
   ## could trap error if required control elements are missing
    nodesToCenter <- control$nodesToCenter
    scale         <- if(!is.null(control$scale))         control$scale         else 1
    ## We could define target in different ways.
    ## I will assume something like target = c("lambda", "gamma")
    ## We could trap errors if the format seems wrong.
    targetNodes <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    numPairs <- length(targetNodes) / 2
    samplers <- nimbleFunctionList(sampler_centered_single_BASE)
    for(i in 1:numPairs) {
       samplers[[i]] <- sampler_centered_single(model, mvSaved, 
                                                c(targetNodes[i], targetNodes[i + numPairs]),
                                                control = list(scale = scale) )
    }
 },
 run = function() {
   mean_covariate <- mean(model[[nodesToCenter]])
   for(i in 1:numPairs) {
      samplers[[i]]$set_mean(mean_covariate)
      samplers[[i]]$run()
   }
 },
 methods = list(
  reset = function () {
   for(i in 1:numPairs) {
     samplers[[i]]$reset()
   }
 })
)


## PAIRED SAMPLER


sampler_RWpaired <- nimbleFunction(
    name = 'sampler_RWpaired',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- if(!is.null(control$log))           control$log           else FALSE
        adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
        adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
        scale         <- if(!is.null(control$scale))         control$scale         else 1
        # reflective    <- if(!is.null(control$reflective))    control$reflective    else FALSE
        ## node list generation
        ## target should be two nodes
        ## first node will get +delta.  second node will get -delta.

        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) != 2)   stop('must use RWpaired sampler on exactly two nodes')
        target1 <- targetAsScalar[1]
        target2 <- targetAsScalar[2]
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory  <- c(0, 0)   ## scaleHistory
        acceptanceHistory  <- c(0, 0)   ## scaleHistory
        if(nimbleOptions('MCMCsaveHistory')) {
            saveMCMChistory <- TRUE
        } else saveMCMChistory <- FALSE
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
##        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(any(model$isDiscrete(targetAsScalar)))
             stop('cannot use RWpaired sampler on discrete-valued target')
##        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
        currentValue1 <- model[[target1]]
        currentValue2 <- model[[target2]]
        propLogScale <- 0
        if(logScale) {
            propLogScale <- rnorm(1, mean = 0, sd = scale)
            propFactor <- exp(propLogScale)
            propValue1 <- currentValue1 * propFactor
            propValue2 <- currentValue2 / propFactor 
        } else {
            shift <- rnorm(1, mean = 0, sd = scale)
            propValue1 <- currentValue1 + shift
            propValue2 <- currentValue2 - shift
        }
        model[[target1]] <<- propValue1
        model[[target2]] <<- propValue2
        logMHR <- calculateDiff(model, target)
        if(logMHR == -Inf) {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            ## Drawing a random number is needed during first testing
            ## of this step in order to keep the random numbers identical
            ## to old behavior to see if tests that depend on particular
            ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
            ## we call decide() to ensure same behavior.
            ## jump <- decide(logMHR)
            ## When new behavior is acceptable, we can remove the above line
            ## and uncomment the following:
            jump <- FALSE
        } else {
            logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf)
            jump <- decide(logMHR)
            if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)         ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale        ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)         ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                }
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                }

                timesRan <<- 0
                timesAccepted <<- 0
        },
        getScaleHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(scaleHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
                return(numeric(1, 0))
            }
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(acceptanceHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
                return(numeric(1, 0))
            }
        },          
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
            }
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)


##notes

assignPairedSampling <- function(mcmcConf, var, scale = 1, removeExistingSamplers = FALSE) {
    nodeNames <- mcmcConf$model$expandNodeNames(var)

  
   if(removeExistingSamplers)    mcmcConf$removeSamplers(var)
    nNodes <- length(nodeNames)
    for(i in 1:(nNodes-1)) {
       mcmcConf$addSampler(type = "RWpaired", target = nodeNames[i:(i+1)], control = list(scale = scale))
    }
}
