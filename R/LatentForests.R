setConstructorS3("LatentForests",
                 function(numLeaves = 0,
                          E = matrix(numeric(0), ncol = 2)) {
                   if (!isBinaryEdgelistVecOrMat(E, numLeaves)) {
                     throw(paste("E is not a valid binary edge list matrix for",
                                 "the given number of leaves ."))
                   }
                   E = edgesInRootedDAGOrder(E)
                   subModels = getSubModelSupports(numLeaves, E)
                   numModels = nrow(subModels)
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   posetAsGraph = subModelsToDAG(subModels)
                   topOrder <- topological.sort(posetAsGraph)
                   tree <- graph.edgelist(E, directed = F)

                   # Parameters
                   dimension = rep(0, numModels)
                   for (j in topOrder) {
                     aE = getSubmodelEdges(subModels[j, ], E)
                     dimension[j] = forestDim(aE, rep(1, nrow(aE)), numLeaves)
                   }

                   extend(
                     ModelPoset(),
                     "LatentForests",
                     .E = E,
                     .subModels = subModels,
                     .numModels = numModels,
                     .prior = prior,
                     .topOrder = topOrder,
                     .tree = tree,
                     .dimension = dimension,
                     .posetAsGraph = posetAsGraph,
                     .numLeaves = numLeaves
                   )
                 })

setMethodS3("getTopOrder", "LatentForests", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

setMethodS3("getPrior", "LatentForests", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

setMethodS3("getNumModels", "LatentForests", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

setMethodS3("getNumLeaves", "LatentForests", function(this) {
  return(this$.numLeaves)
}, appendVarArgs = F)

setMethodS3("getNumSamples", "LatentForests", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

setMethodS3("getNumVertices", "LatentForests", function(this) {
  numLeaves = this$getNumLeaves()
  if (numLeaves == 0) {
    return(0)
  }
  return(max(2 * this$getNumLeaves() - 2, 1))
}, appendVarArgs = F)

setMethodS3("getModelWithSupport", "LatentForests", function(this, support) {
  if (length(support) != nrow(this$getAllEdges())) {
    throw("Invalid support length.")
  }
  if (is.null(this$.subModelsAsStrings)) {
    this$.subModelsAsStrings = apply(this$.subModels, 1, function(x) { paste(x, collapse="") })
  }
  return(which(this$.subModelsAsStrings == paste(support, collapse = "")))
}, appendVarArgs = F)

setMethodS3("getDimension", "LatentForests", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)

setMethodS3("setData", "LatentForests", function(this, X) {
  this$.X = X
  this$.sampleCovMat = t(X) %*% X
}, appendVarArgs = F)

setMethodS3("getData", "LatentForests", function(this) {
  if (is.null(this$.X)) {
    throw("No data has been set for models.")
  }
  return(this$.X)
}, appendVarArgs = F)

setMethodS3("getSamplingCovMat", "LatentForests", function(this) {
  if (is.null(this$.X)) {
    throw("No data has been set for models.")
  }
  return(this$.sampleCovMat)
}, appendVarArgs = F)

setMethodS3("parents", "LatentForests", function(this, model) {
  if (length(model) != 1) {
    throw("parents can only accept a single model.")
  }
  return(as.numeric(neighbors(this$posetAsGraph, model, "in")))
}, appendVarArgs = F)

setMethodS3("logLikeMle", "LatentForests", function(this, model) {
  return((this$emMain(model, starts = 5, maxIter = 1000, tol = 1e-4))$logLike)
}, appendVarArgs = F)

setMethodS3("getSupport", "LatentForests", function(this, model) {
  return(this$.subModels[model, ])
}, appendVarArgs = F)

setMethodS3("getAllEdges", "LatentForests", function(this, model) {
  return(this$.E)
}, appendVarArgs = F)

setMethodS3("logLike", "LatentForests", function(this, covMat) {
  n = this$getNumSamples()
  tXX = this$getSamplingCovMat()
  return(
    - n / 2 * ncol(covMat) * log(2 * pi)
    - n / 2 * as.numeric(determinant(covMat)$modulus)
    - (1 / 2) * sum(tXX * chol2inv(chol(covMat)))
  )
}, appendVarArgs = F)

setMethodS3("emMain", "LatentForests", function(this, model, starts = 5,
                                                maxIter = 1000, tol = 1e-4) {
  bestLogLike <- -Inf
  n = this$getNumSamples()
  numVertices = this$getNumVertices()
  support = this$getSupport(model)
  numLeaves = this$getNumLeaves()
  if (numLeaves == 1) {
    return(list(logLike = this$logLike(matrix(1)), covMat = matrix(1)))
  } else if (numLeaves == 2) {
    if (sum(support) == 0) {
      return(list(logLike = this$logLike(diag(2)), covMat = diag(2)))
    } else {
      covMat = cor(this$getData())
      return(list(logLike = this$logLike(covMat), covMat = covMat))
    }
  }

  for (run in 1:starts) {
    ## Initializing correlations to be Uniform(.2,.9)
    edgeCorrelations = (runif(length(support), 0.2, .9)) * support
    curCovMat = this$getCovMat(edgeCorrelations)

    curLogLike = this$logLike(curCovMat[1:numLeaves, 1:numLeaves, drop = F])

    conv = FALSE
    iter = 1
    while (conv == FALSE & iter < maxIter) {
      curCovMat = this$emSteps(support, curCovMat)
      if (iter < 10 || iter %% 5 == 0) {
        nextLogLike = this$logLike(curCovMat[1:numLeaves, 1:numLeaves, drop = F])
        if (1 - nextLogLike / curLogLike < tol) {
          conv = TRUE
        }
        curLogLike = nextLogLike
      }
      iter = iter + 1
    }

    if (iter >= maxIter) {
      print("Maximum EM Iterations used!")
    }
    if (curLogLike > bestLogLike) {
      bestLogLike = curLogLike
      bestCovMat = curCovMat
    }
  }
  return(list(logLike = bestLogLike, covMat = bestCovMat))
}, appendVarArgs = F)

setMethodS3("getCovMat", "LatentForests", function(this, edgeCorrelations) {
  # for given edge correlations it output the full correlation matrix
  # (!!) this works only if edges were provided so that
  # T forms a rooted tree (we use the DAG parametrization)
  v = this$getNumVertices()
  E = this$getAllEdges()

  if (nrow(E) == 0) {
    return(diag(v))
  }

  Lam = matrix(0, v, v)
  U = diag(v)

  for (i in 1:length(edgeCorrelations)) {
    e1 = E[i, 1]
    e2 = E[i, 2]
    Lam[e1, e2] = edgeCorrelations[i]
    U[e2, e2] = 1 - edgeCorrelations[i] ^ 2
  }
  M = solve(diag(v) - Lam)
  return(t(M) %*% U %*% M)
}, appendVarArgs = F)

setMethodS3("emSteps", "LatentForests", function(this, support, S) {
  tXX = this$getSamplingCovMat()
  E = this$getAllEdges()
  v = nrow(S) # Num vertices
  m = nrow(tXX) # Num leaves
  n = this$getNumSamples()

  SmmInvtXX = solve(S[1:m, 1:m, drop = F], tXX)
  Sigma11 = 1 / n * tXX
  Sigma12 = 1 / n * t(S[(m + 1):v, 1:m, drop = F] %*% SmmInvtXX)

  SmmInvSmv = solve(S[1:m, 1:m, drop = F], S[1:m, (m + 1):v, drop = F])
  schurComplement = S[(m + 1):v, (m + 1):v, drop = F] - S[(m + 1):v, 1:m, drop = F] %*% SmmInvSmv
  Sigma22 = schurComplement + 1 / n * S[(m + 1):v, 1:m, drop = F] %*% SmmInvtXX %*% SmmInvSmv
  Sigma = rbind(cbind(Sigma11, Sigma12), cbind(t(Sigma12), Sigma22))

  edgeCorrelations = rep(0, length(support))
  if (nrow(E) > 0) {
    for (i in 1:(nrow(E))) {
      if (support[i] == 1) {
        e1 = E[i, 1]
        e2 = E[i, 2]
        edgeCorrelations[i] = Sigma[e1, e2] / sqrt(Sigma[e1, e1] * Sigma[e2, e2])
      }
    }
  }
  return(this$getCovMat(edgeCorrelations))
}, appendVarArgs = F)

# This function must compute the learning coefficient of subModel in superModel
# Returns a named list with components lambda and m
setMethodS3("learnCoef", "LatentForests", function(this, superModel, subModel) {
  support = this$getSupport(superModel)
  subSupport = this$getSupport(subModel)
  E = this$getAllEdges()
  numLeaves = this$getNumLeaves()

  w.sum = 0
  m = 1
  if (!all(support == subSupport)) {
    # Compute w sum (primary lambda component)
    subEdges = E[subSupport == 1, , drop = F]
    qForestNodes = union(subEdges, 1:numLeaves)

    for (it in 1:length(support)) {
      if (support[it] - subSupport[it] == 1) {
        w.sum = w.sum + sum(E[it, ] %in% qForestNodes)
      }
    }

    # Compute multiplicity
    nodeDegsBig = table(E[support == 1, , drop = F])
    deg2NodesBig = as.numeric(names(nodeDegsBig[nodeDegsBig == 2]))
    m = 1 + length(setdiff(deg2NodesBig, qForestNodes))
  }

  lambda = this$getDimension(subModel) / 2 + w.sum / 4

  return(list(lambda = lambda, m = m))
}, appendVarArgs = F)
