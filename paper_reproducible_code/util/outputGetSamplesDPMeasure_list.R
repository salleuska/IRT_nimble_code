rm(list=ls())
library(nimble)

## FUNCTION THAT WE NEED
setupResults <- function() { # Use lexical scoping
  myResults <- list()
  setMatrix <- function(m, i) {
    myResults[[i]] <<- m
  }
  getResults <- function() {
    myResults
  }
  ## We need two things in the R global environment.
  ## I will give them hidden names.
  ## The nimbleRcall system could be more general for scoping, but it currently isn't.
  .GlobalEnv$.setMatrix_for_samplesDPmeasure_internal <- setMatrix
  .GlobalEnv$.setMatrix_for_samplesDPmeasure <- nimbleRcall(prototype = function(m = double(2), i = double()) {},
                                                            returnType = void(),
                                                            Rfun = '.setMatrix_for_samplesDPmeasure_internal')
  list(getResults = getResults)
}

getDPmeasureResults <- setupResults()

findClusterNodes <- function(model, target) {
  ## Determine which model nodes are the cluster parameters by processing expressions to look
  ## for what is indexed by the dCRP clusterID nodes. This also determine which clusterID
  ## each cluster parameter is associated with.
  targetVar <- model$getVarNames(nodes = target)
  targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  deps <- model$getDependencies(target, self = FALSE)
  declIDs <- sapply(deps, function(x) model$getDeclID(x))
  uniqueIDs <- unique(declIDs)
  depsByDecl <- lapply(uniqueIDs, function(x) deps[which(x == declIDs)])
  ## Find one example dependency per BUGS declaration for more efficient processing
  exampleDeps <- sapply(depsByDecl, `[`, 1)
  
  ## Once we find the cluster parameter variables below, we want to evaluate the cluster membership
  ## values (e.g., xi[1],...,xi[n]) for all possible values they could take, this will
  ## allow us to determine all possible cluster nodes in the model (though some may
  ## not actually be specified in the model, if there is truncation).
  ## Therefore, set up an evaluation environment in which (xi[1],...,xi[n]) = (1,2,...,n)
  ## first try was: e[[targetVar]] <- seq_along(targetElements)
  ## However in first try, that wouldn't handle xi[3:10] ~ dCRP(), but next construction does.
  e <- list()
  idxExpr <- model$getDeclInfo(target)[[1]]$indexExpr[[1]]
  eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(targetElements)), list(VAR = targetVar, IDX = idxExpr)))
  ## For cases of cross clustering (e.g., mu[xi[i],eta[j]]) we need the other dcrp node(s)
  nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  dists <- model$getDistribution(nodes)
  if(length(dists == 'dCRP') > 1) { 
    dcrpNodes <- nodes[dists == 'dCRP' & nodes != target]
    for(i in seq_along(dcrpNodes)) {
      dcrpElements <- model$expandNodeNames(dcrpNodes[i], returnScalarComponents = TRUE)
      dcrpVar <- model$getVarNames(nodes = dcrpNodes[i])
      idxExpr <- model$getDeclInfo(dcrpNodes[i])[[1]]$indexExpr[[1]]
      eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(dcrpElements)), list(VAR = dcrpVar, IDX = idxExpr)))
      
    }
  }
  
  clusterNodes <- indexExpr <- clusterIDs <- list()
  clusterVars <- indexPosition <- numIndexes <- targetIsIndex <- targetIndexedByFunction <-
    loopIndex <- NULL
  varIdx <- 0
  
  targetNonIndex <- NULL
  multipleStochIndexes <- NULL
  
  modelVars <- model$getVarNames()
  modelVars <- modelVars[!modelVars == targetVar]
  
  ## Process model declaration expressions to find stochastic indexing and the indexed variable.
  for(idx in seq_along(exampleDeps)) {
    ## Pull out expressions, either as RHS of deterministic or parameters of stochastic
    fullExpr <- cc_getNodesInExpr(model$getValueExpr(exampleDeps[idx]))
    for(j in seq_along(fullExpr)) {
      subExpr <- parse(text = fullExpr[j])[[1]]  # individual parameter of stochastic or RHS of deterministic
      len <- length(subExpr)
      ## Look for target variable within expression, but only when used within index
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' &&
         sum(all.vars(subExpr) == targetVar) && subExpr[[2]] != targetVar) {
        varIdx <- varIdx + 1
        multipleStochIndexes <- c(multipleStochIndexes, FALSE)
        
        clusterVars <- c(clusterVars, deparse(subExpr[[2]]))
        
        ## Determine which index the target variable occurs in.
        k <- whichIndex <- 3
        foundTarget <- FALSE
        while(k <= len) {
          if(sum(all.vars(subExpr[[k]]) == targetVar)) {
            if(foundTarget) {
              stop("findClusterNodes: CRP variable used multiple times in ", deparse(subExpr),
                   ". NIMBLE's CRP MCMC sampling not designed for this situation.")
            } else {
              foundTarget <- TRUE
              whichIndex <- k
            }
          }
          ## We will need to relax this when allow crossed clustering.
          if(sum(all.vars(subExpr[[k]]) %in% modelVars)) { ## cases like mu[xi[i],eta[j]]
            ## We are adding support for this case.
            ## warning("findClusterNodes: multiple indexing variables in '", deparse(subExpr),
            ##          "'. NIMBLE's CRP MCMC sampling not designed for this situation.")
            multipleStochIndexes[varIdx] <- TRUE
          }                
          k <- k+1
        }
        if(!foundTarget) stop("findClusterNodes: conflicting information about presence of CRP variable in expression.")
        
        declInfo <-  model$getDeclInfo(exampleDeps[idx])[[1]]
        
        ## Determine how target variable enters into cluster node definition
        indexPosition[varIdx] <- whichIndex-2
        numIndexes[varIdx] <- len - 2
        indexExpr[[varIdx]] <- subExpr
        ## Is target used directly as index, e.g., "mu[xi[.]]" as opposed to something like "mu[1+xi[.]]".
        targetIsIndex[varIdx] <- length(subExpr[[whichIndex]]) == 3 &&
          subExpr[[whichIndex]][[1]] == '[' &&
          subExpr[[whichIndex]][[2]] == targetVar
        ## Is indexing of target a simple index, e.g. xi[i], as opposed to something like "xi[n-i+1]".
        targetIndexedByFunction[varIdx] <- any(sapply(declInfo$symbolicParentNodes,
                                                      function(x) 
                                                        length(x) >= 3 && x[[1]] == '[' &&
                                                        x[[2]] == targetVar && length(x[[3]]) > 1))
        ## Determine all sets of index values so they can be evaluated in context of possible values of target element values.
        unrolledIndices <- declInfo$unrolledIndicesMatrix
        
        if(targetIndexedByFunction[varIdx] && ncol(unrolledIndices) > 1)  ## Now that we allow cluster parameters with multiple indexes, this is very hard to handle in terms of identifying what column of unrolledIndices to use for sorting clusterNodes.
          stop("findClusterNodes: Detected that a cluster parameter is indexed by a function such as 'mu[xi[n-i+1]]' rather than simple indexing such as 'mu[xi[i]]'. NIMBLE's CRP MCMC sampling not designed for this case.")
        loopIndexes <- unlist(sapply(declInfo$symbolicParentNodes,
                                     function(x) {
                                       if(length(x) >= 3 && x[[1]] == '[' &&
                                          x[[2]] == targetVar) return(deparse(x[[3]]))
                                       else return(NULL) }))
        if(length(loopIndexes) != 1)
          stop("findClusterNodes: found cluster membership parameters that use different indexing variables; NIMBLE's CRP sampling not designed for this case.")
        ## Note not clear when NULL would be the result...
        loopIndex[[varIdx]] <- loopIndexes
        
        ## Determine potential cluster nodes by substituting all possible clusterID values into the indexing expression. 
        n <- nrow(unrolledIndices)
        if(n > 0 && loopIndex %in% dimnames(unrolledIndices)[[2]]) {  # catch cases like use of xi[2] rather than xi[i]
          ## Order so that loop over index of cluster ID in order of cluster ID so that
          ## clusterNodes will be grouped in chunks of unique cluster IDs for correct
          ## sampling of new clusters when have multiple obs per cluster.
          ord <- order(unrolledIndices[ , loopIndex[varIdx]])
          unrolledIndices <- unrolledIndices[ord, , drop = FALSE]
          clusterIDs[[varIdx]] <- unrolledIndices[ , loopIndex[varIdx]]
          
          clusterNodes[[varIdx]] <- rep(NA, n)
          
          ## Determine unevaluated expression, e.g., muTilde[xi[i],j] not muTilde[xi[1],2]
          expr <- declInfo$valueExprReplaced
          expr <- parse(text = cc_getNodesInExpr(expr)[[j]])[[1]]
          templateExpr <- expr   
          
          ## Now evaluate index values for all possible target element values, e.g.,
          ## xi[i] for all 'i' values with xi taking values 1,...,n
          for(i in seq_len(n)) { 
            for(k in 3:len) # this will deal with muTilde[xi[i], j] type cases
              if(length(all.vars(expr[[k]])))  # prevents eval of things like 1:3, which the as.numeric would change to c(1,3)
                templateExpr[[k]] <- as.numeric(eval(substitute(EXPR, list(EXPR = expr[[k]])),
                                                     c(as.list(unrolledIndices[i,]), e)))  # as.numeric avoids 1L, 2L, etc.
              clusterNodes[[varIdx]][i] <- deparse(templateExpr)  # convert to node names
          }
        } else {
          clusterNodes[[varIdx]] <- character(0)
          clusterIDs[[varIdx]] <- numeric(0)
        }
      } 
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' && subExpr[[2]] == targetVar)
        targetNonIndex <- deparse(model$getDeclInfo(exampleDeps[idx])[[1]]$codeReplaced)
    }
  }
  
  ## Find the potential cluster nodes that are actually model nodes,
  ## making sure that what we decide are real cluster nodes are the full potential set
  ## or a truncated set that starts with the first cluster node, e.g., muTilde[1], ..., muTilde[3] is ok;
  ## muTilde[2], ..., muTilde[4] is not (unless the model nodes are muTilde[2], ...., muTilde[4]).
  nTilde <- sapply(clusterNodes, length)
  modelNodes <- model$getNodeNames()
  
  for(varIdx in seq_along(clusterVars)) {
    if(nTilde[varIdx]) {
      if(any(is.na(clusterNodes[[varIdx]])))  
        stop("findClusterNodes: fewer cluster IDs in ", target, " than elements being clustered.")
      
      ## Handle cases where indexing of variables in dynamic indexing does not correspond to actual
      ## stochastic model nodes.
      if(any(!clusterNodes[[varIdx]] %in% modelNodes)) {
        tmp <- mapply(function(node, id) {
          if(!node %in% modelNodes) {
            node <- model$expandNodeNames(node)
            id <- rep(id, length(node))
          }
          return(list(nodes = node, ids = id))
        }, clusterNodes[[varIdx]], clusterIDs[[varIdx]])
        dimnames(tmp) <- NULL
        clusterNodes[[varIdx]] <- unlist(tmp[1, ])
        clusterIDs[[varIdx]] <- unlist(tmp[2, ])
      }
      ## Now remove duplicates when indexed variables correspond to same model node,
      ## but only for duplicates within a cluster.
      groups <- split(clusterNodes[[varIdx]], clusterIDs[[varIdx]])
      dups <- unlist(lapply(groups, duplicated))
      clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][!dups]
      clusterIDs[[varIdx]] <- clusterIDs[[varIdx]][!dups]
      
      ## Formerly we were checking that we had a contiguous set of cluster nodes
      ## starting with the first one, but for clusterNodes with more than one index and
      ## truncation this is hard to do, so just fall back to returning the clusterNodes
      ## that are actually part of the model.
      validNodes <- clusterNodes[[varIdx]] %in% modelNodes
      
      if(!all(validNodes)) {  # i.e., truncated representation
        clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][validNodes]
        clusterIDs[[varIdx]] <- clusterIDs[[varIdx]][validNodes]
      }
    }
  }
  
  
  nTilde <- sapply(clusterNodes, length)
  numNodesPerCluster <- sapply(clusterIDs, function(x) {
    tbl <- table(x)
    num <- unique(tbl)
    if(length(num) > 1) stop("findClusterNodes: detected differing numbers of nodes (i.e., parameters) per cluster. NIMBLE's CRP sampling not designed for this case.")
    return(num)})
  
  return(list(clusterNodes = clusterNodes, clusterVars = clusterVars, nTilde = nTilde,
              numNodesPerCluster = numNodesPerCluster, clusterIDs = clusterIDs, loopIndex = loopIndex,
              targetIsIndex = targetIsIndex, indexPosition = indexPosition, indexExpr = indexExpr,
              numIndexes = numIndexes, targetIndexedByFunction = targetIndexedByFunction,
              targetNonIndex = targetNonIndex, multipleStochIndexes = multipleStochIndexes))
}

getSamplesDPmeasureNames <- function(clusterVarInfo, model, truncG, p) {
  result <- NULL
  for(j in 1:p) {
    tildeNodesModel <- model$expandNodeNames(clusterVarInfo$clusterVars[j], returnScalarComponents=TRUE) # tilde nodes j in model
    allIndexes <- 1:length(tildeNodesModel)
    
    clusterID <- 1
    tildeNodesPerClusterID <- model$expandNodeNames(clusterVarInfo$clusterNodes[[j]][clusterVarInfo$clusterIDs[[j]] == clusterID], returnScalarComponents=TRUE) # tilde nodes in cluster with id 1
    aux <- match(tildeNodesModel, tildeNodesPerClusterID, nomatch = 0) 
    cn <-  tildeNodesModel[which(aux != 0)]
    
    tmp <- matrix('', length(cn), truncG)
    for(i in seq_along(cn)) {   # for each element of the cluster parameters, replicate 1:truncG times
      expr <- parse(text = cn[i])[[1]]
      tmp[i, ] <- sapply(seq_len(truncG),
                         function(idx) {
                           expr[[2+clusterVarInfo$indexPosition[j]]] <- as.numeric(idx)
                           return(deparse(expr))
                         })
    }
    result <- rbind(result , tmp)
  }
  return(c(result))
}


# thefollowing is code we currently use
getSamplesDPmeasure_list <- function(MCMC, epsilon = 1e-4) {
  if(exists('model',MCMC, inherits = FALSE)) compiled <- FALSE else compiled <- TRUE
  if(compiled) {
    if(!exists('Robject', MCMC, inherits = FALSE) || !exists('model', MCMC$Robject, inherits = FALSE))
      stop("getSamplesDPmeasure: problem with finding model object in compiled MCMC")
    model <- MCMC$Robject$model
    mvSamples <- MCMC$Robject$mvSamples
  } else {
    model <- MCMC$model
    mvSamples <- MCMC$mvSamples
  }
  
  getSamplesDPmeasureC <- sampleDPmeasure_list(model, mvSamples, epsilon)
  
  if(compiled) {
    CgetSamplesDPmeasureC <- compileNimble(getSamplesDPmeasureC, project = model)
    CgetSamplesDPmeasureC$run()
  } else {
    getSamplesDPmeasureC$run()
  }
  samplesMeasure <- getDPmeasureResults$getResults()
  
  dcrpVar <- getSamplesDPmeasureC$dcrpVar
  clusterVarInfo <- findClusterNodes(model, dcrpVar) 
  namesVars <- getSamplesDPmeasureC$tildeVars
  p <- length(namesVars)
  namesVars <- getSamplesDPmeasureNames(clusterVarInfo, model, 1, p)
  
  niter <- length(samplesMeasure)
  
  for(i in 1:niter) {
    colnames(samplesMeasure[[i]]) <- c("weights", namesVars)
  }
  
  return(samplesMeasure)
}




sampleDPmeasure_list <- nimbleFunction(
  name = 'sampleDPmeasure_list',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, epsilon){
    mvSavedVars <- mvSaved$varNames
    
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    distributions <- model$getDistribution(stochNodes) 
    
    ## Determine if there is a dCRP-distributed node and that it is monitored.
    dcrpIndex <- which(distributions == 'dCRP')
    if(length(dcrpIndex) == 1) {
      dcrpNode <- stochNodes[dcrpIndex] 
      dcrpVar <- model$getVarNames(nodes = dcrpNode)
    } else {
      if(length(dcrpIndex) == 0 ){
        stop('sampleDPmeasure: One node with a dCRP distribution is required.\n')
      }
      stop('sampleDPmeasure: Currently only models with one node with a dCRP distribution are allowed.\n')
    }
    if(sum(dcrpVar == mvSavedVars) == 0)
      stop('sampleDPmeasure: The node having the dCRP distribution has to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    
    ## Find the cluster variables, named tildeVars
    dcrpElements <- model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE)
    clusterVarInfo <- nimble:::findClusterNodes(model, dcrpVar) 
    tildeVars <- clusterVarInfo$clusterVars
    if( is.null(tildeVars) )  ## probably unnecessary as checked in CRP sampler, but best to be safe
      stop('sampleDPmeasure: The model should have at least one cluster variable.\n')
    
    ## Check that cluster parameters are IID (across clusters, non-IID within clusters is ok), as required for random measure G
    isIID <- TRUE
    for(i in seq_along(clusterVarInfo$clusterNodes)) {
      clusterNodes <- clusterVarInfo$clusterNodes[[i]]  # e.g., 'thetatilde[1]',...,
      clusterIDs <- clusterVarInfo$clusterIDs[[i]]
      splitNodes <- split(clusterNodes, clusterIDs)
      valueExprs <- lapply(splitNodes, function(x) {
        out <- sapply(x, model$getValueExpr)
        names(out) <- NULL
        out
      })
      if(length(unique(valueExprs)) != 1) 
        isIID <- FALSE
    }
    
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvGammaConjugacy(model, clusterVarInfo, length(dcrpElements), 'dinvgamma'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvGammaConjugacy(model, clusterVarInfo, length(dcrpElements), 'dgamma'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvWishartConjugacy(model, clusterVarInfo, length(dcrpElements), 'dinvwish'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvWishartConjugacy(model, clusterVarInfo, length(dcrpElements), 'dwish'))
      isIID <- TRUE
    ## Tricky as MCMC might not be using conjugacy, but presumably ok to proceed regardless of how
    ## MCMC was done, since conjugacy existing would guarantee IID.
    if(!isIID) stop('sampleDPmeasure: cluster parameters have to be independent and identically distributed. \n')
    
    ## Check that necessary variables are being monitored.
    
    ## Check that cluster variables are monitored.
    counts <- tildeVars %in% mvSavedVars
    if( sum(counts) != length(tildeVars) ) 
      stop('sampleDPmeasure: The node(s) representing the cluster variables must be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
    
    parentNodesTildeVars <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes %in% unlist(clusterVarInfo$clusterNodes)]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE)
      for(j in seq_along(tildeVars)) {
        if(sum(aux == clusterVarInfo$clusterNodes[[j]][1]))
          parentNodesTildeVars <- c(parentNodesTildeVars, candidateParentNodes[i])
      }
    }
    if(length(parentNodesTildeVars)) {
      parentNodesTildeVarsDeps <- model$getDependencies(parentNodesTildeVars, self = FALSE)
    } else parentNodesTildeVarsDeps <- NULL
    ## make sure tilde nodes are included (e.g., if a tilde node has no stoch parents) so they get simulated
    parentNodesTildeVarsDeps <- model$topologicallySortNodes(c(parentNodesTildeVarsDeps, unlist(clusterVarInfo$clusterNodes)))
    
    if(!all(model$getVarNames(nodes = parentNodesTildeVars) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the cluster variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    if(is.null(parentNodesTildeVars)) parentNodesTildeVars <- tildeVars  ## to avoid NULL which causes compilation issues
    
    ## Check that parent nodes of cluster IDs are monitored.   
    parentNodesXi <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE)
      if(sum(aux == dcrpNode)) {
        parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
      }
    }
    
    if(!all(model$getVarNames(nodes = parentNodesXi) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the membership variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    if(is.null(parentNodesXi)) parentNodesXi <- dcrpNode  ## to avoid NULL which causes compilation issues
    
    ## End of checks of monitors.
    
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    if(length(parentNodesXi)) {
      fixedConc <- FALSE
      parentNodesXiDeps <- model$getDependencies(parentNodesXi, self = FALSE)
      parentNodesXiDeps <- parentNodesXiDeps[!parentNodesXiDeps == dcrpNode]
    } else {
      parentNodesXiDeps <- dcrpNode
    }
    
    dataNodes <- model$getDependencies(dcrpNode, stochOnly = TRUE, self = FALSE)
    N <- length(model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE))
    
    p <- length(tildeVars)
    lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
    dimTildeVarsNim <- numeric(p+1) # nimble dimension (0 is scalar, 1 is 2D array, 2 is 3D array) (dimTildeVarsNim=dimTildeNim)
    dimTildeVars <- numeric(p+1) # dimension to be used in run code (dimTildeVars=dimTilde)
    for(i in 1:p) {
      dimTildeVarsNim[i] <- model$getDimension(clusterVarInfo$clusterNodes[[i]][1])
      dimTildeVars[i] <- lengthData^(dimTildeVarsNim[i]) 
    }
    nTildeVarsPerCluster <-  clusterVarInfo$numNodesPerCluster
    nTilde <- numeric(p+1)
    nTilde[1:p] <- clusterVarInfo$nTilde / nTildeVarsPerCluster
    if(any(nTilde[1:p] != nTilde[1])){
      stop('sampleDPmeasure: All cluster parameters must have the same number of parameters.\n')
    }
    
    
    tildeVarsCols <- c(dimTildeVars[1:p]*nTildeVarsPerCluster, 0)
    tildeVarsColsSum <- c(0, cumsum(tildeVarsCols))
    mvIndexes <- matrix(0, nrow=nTilde[1], ncol=(sum(dimTildeVars[1:p]*nTildeVarsPerCluster))) 
    for(j in 1:p) {
      tildeNodesModel <- model$expandNodeNames(clusterVarInfo$clusterVars[j], returnScalarComponents=TRUE) # tilde nodes j in model
      allIndexes <- 1:length(tildeNodesModel)
      for(l in 1:nTilde[1]) {
        clusterID <- l
        tildeNodesPerClusterID <- model$expandNodeNames(clusterVarInfo$clusterNodes[[j]][clusterVarInfo$clusterIDs[[j]] == clusterID], returnScalarComponents=TRUE) # tilde nodes in cluster with id 1
        aux <- match(tildeNodesModel, tildeNodesPerClusterID, nomatch = 0) 
        mvIndexes[l,(tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1] ] <- which(aux != 0)
      }
    }
    
    
    ## Storage object to be sized in run code based on MCMC output.
    #samples <- matrix(0, nrow = 1, ncol = 1)   
    ## Truncation level of the random measure 
    #truncG <- 0 
    
    ## control list extraction
    ## The error of approximation G is given by (conc / (conc +1))^{truncG-1}. 
    ## we are going to define an error of approximation and based on the posterior values of the conc parameter define the truncation level of G
    ## the error is between errors that are considered very very small in the folowing papers
    ## Ishwaran, H., & James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association, 96(453), 161-173.
    ## Ishwaran, H., & Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. Biometrika, 87(2), 371-390.
    # epsilon <- 1e-4
    niter <- 0
    setupOutputs(lengthData, dcrpVar)
    
  },
  
  run=function(){
    
    niter <<- getsize(mvSaved) # number of iterations in the MCMC
    
    # getting posterior samples of concentartion parameter
    if( fixedConc ) { # we don't need the output from the mcmc
      concSamples <- nimNumeric(length = niter, value = model$getParam(dcrpNode, 'conc'))
    } else { # we need the output from the mcmc
      concSamples <- numeric(niter)
      for( iiter in 1:niter ) {
        nimCopy(from = mvSaved, to = model, nodes = parentNodesXi, row=iiter) 
        model$calculate(parentNodesXiDeps)
        concSamples[iiter] <- model$getParam(dcrpNode, 'conc')
      }
    }
    
    # posterior samples of random measure
    for(iiter in 1:niter){
      checkInterrupt()
      
      # truncation level (to be reduced later)
      concAux <- concSamples[iiter] + N
      truncG <- log(epsilon) / log(concAux / (concAux+1)) 
      truncG <- ceiling(truncG) 
      
      ## sampling stick-breaking weights: as many as truncG
      weigthsAux <- nimNumeric(truncG)
      vaux <- rbeta(1, 1, concAux)
      v1prod <- 1
      weigthsAux[1] <- vaux  
      for(l1 in 2:truncG) {
        v1prod <- v1prod * (1-vaux)
        vaux <- rbeta(1, 1, concAux)
        weigthsAux[l1] <- vaux * v1prod 
      }
      weigthsAux[1:truncG] <- weigthsAux[1:truncG] / (1 - v1prod * (1-vaux)) # normalizing weigths
      
      ## getting the unique cluster variables that where sampled in the mcmc and the sampling probabilities of the polya urn of the unique cluster variables
      probs <- nimNumeric(N) # polya urn probabilities
      uniqueValues <- matrix(0, nrow = N, ncol = tildeVarsColsSum[p+1])  # unique cluster variables
      xiiter <- mvSaved[dcrpVar, iiter]
      range <- min(xiiter):max(xiiter) 
      index <- 1
      for(i in seq_along(range)){   
        cond <- sum(xiiter == range[i])
        if(cond > 0){
          probs[index] <- cond
      ## slight workaround because can't compile mvSaved[tildeVars[j], iiter]  
          nimCopy(mvSaved, model, tildeVars, row = iiter)
          for(j in 1:p){
            jcols <- (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]
            uniqueValues[index, jcols] <- values(model, tildeVars[j])[mvIndexes[range[i], jcols]]
          }
          index <- index+1
        }
      }
      probs[index] <- concSamples[iiter] # last probability of the polya urn (the sample concetration parameter)
      newValueIndex <- index 

      
      ## copy tilde parents into model for use in simulation below when simulate atoms of G_0  
      nimCopy(mvSaved, model, parentNodesTildeVars, row = iiter)
      
      ## reduced samples from random measure: 
      indexesG <- nimNumeric(truncG) # indicates if an existing or a new atom is sampled from the polya urn. New tom are indicated by newValueIndex 
      weigths <- nimNumeric(truncG) # weights if random measure
      atoms <- matrix(0, ncol = tildeVarsColsSum[p+1], nrow = truncG) # atoms if random measure
      indexG0 <- newValueIndex # used for new atoms
      for(l1 in 1:truncG) {
        indexesG[l1] <- rcat(prob = probs[1:newValueIndex])
        if(indexesG[l1] < newValueIndex) { # an existing atom was sampled and the corresponding weights are added
          weigths[indexesG[l1]] <- weigths[indexesG[l1]] + weigthsAux[l1] 
          sumCol <- 0
          for(j in 1:p){
            jcols <- (sumCol + 1):(sumCol + tildeVarsCols[j]) 
            atoms[indexesG[l1], jcols] <-  uniqueValues[indexesG[l1], (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]]  
            sumCol <- sumCol + tildeVarsCols[j] 
          }
        } else { # a new atom is sampled (from G0) with its corresponding new weight
          weigths[indexG0] <-  weigthsAux[l1] 
          model$simulate(parentNodesTildeVarsDeps)
          sumCol <- 0
          for(j in 1:p){
            jcols <- (sumCol + 1):(sumCol + tildeVarsCols[j])
            atoms[indexG0, jcols] <- values(model, tildeVars[j])[mvIndexes[1,  (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]]] # <<-
            sumCol <- sumCol + tildeVarsCols[j] 
          }
          indexG0 <- indexG0 + 1
        }
      }
      # the random measure has indexG0-1 atoms
      
      # check that all unique tilde variables were actually sampled. If not we need to rearrange when creating final output
      missingIndex <- nimNumeric(newValueIndex-1)
      uniqueIndex <- 1:(newValueIndex - 1)
      ii <- 1
      for(i in seq_along(uniqueIndex)) {
        if( !any(indexesG == uniqueIndex[i]) ) {
          missingIndex[ii] <- uniqueIndex[i]
          ii <- ii + 1
        }
      }
      
      # final output
      nrow0 <- indexG0 - 1 + ii - 1
      output <- nimNumeric(nrow0*(tildeVarsColsSum[p+1]+1))
      ii <- 1
      imissing <- 1 
      for(i in 1:(indexG0 - 1)) { # saving weights
        if(i != missingIndex[imissing]) {
          output[ii] <- weigths[i]
          ii <- ii + 1  
        } else {
          imissing <- imissing + 1
        }
      }
      for(j in 1:tildeVarsColsSum[p+1]) { # saving atoms
        imissing <- 1
        for(i in 1:(indexG0 - 1)) {
          if(i != missingIndex[imissing]) {
            output[ii] <- atoms[i, j]
            ii <- ii + 1  
          } else {
            imissing <- imissing + 1
          }
        }
      }
      
      
      output <- nimNumeric(value = c(weigths[1:nrow0], atoms[1:nrow0, 1:tildeVarsColsSum[p+1]]), length=nrow0*(tildeVarsColsSum[p+1]+1))
      
      sampleG = matrix(value = output, nrow = nrow0, ncol = (tildeVarsColsSum[p+1]+1), type = "double") #(p+1)
      .setMatrix_for_samplesDPmeasure(sampleG, iiter)
    }
  },
  
  methods = list( reset = function () {} )
)

#-----------------------------------------
# EXAMPLE 
#-----------------------------------------
code <- nimbleCode({
  xi[1:n] ~ dCRP(alpha, n)
  for(i in 1:n){
    mu[i] ~ dnorm(0, var=100)
    y[i] ~ dnorm(mu[xi[i]],  var=1)
  }
  alpha ~ dgamma(1, 1)
})
n <- 100
Consts <- list(n = n)
set.seed(1)
Inits <- list(xi = sample(1:10, size=n, replace=TRUE), 
              mu = rnorm(n, 0, sd=sqrt(100)),
              alpha = 1) 
Data <- list(y = c(rnorm(n/2, -5, 1), rnorm(n/2, 5, 1)))
m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
cm <- compileNimble(m)
mConf <- configureMCMC(m, monitors = c('xi','mu', 'alpha'), print=TRUE)  
mMCMC <- buildMCMC(mConf)
cMCMC <- compileNimble(mMCMC, project = m)
cMCMC$run(100)

aux <- getSamplesDPmeasure_list(cMCMC, epsilon = 1e-4)
aux



