setConstructorS3("FactorAnalyses",
                 function(numCovariates = 1, maxNumFactors = 0) {
                   numModels = maxNumFactors + 1
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (maxNumFactors == 0) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(1, numModels)
                   for (j in topOrder) {
                     k = j - 1
                     dimension[j] = (k + 1) * numCovariates - choose(k, 2)
                   }

                   extend(
                     ModelPoset(),
                     "FactorAnalyses",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .numCovariates = numCovariates,
                     .maxNumFactors = maxNumFactors
                   )
                 })

setMethodS3("getTopOrder", "FactorAnalyses", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

setMethodS3("getPrior", "FactorAnalyses", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

setMethodS3("getNumModels", "FactorAnalyses", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

setMethodS3("getNumFactorsForModel", "FactorAnalyses", function(this, model) {
  return(model - 1)
}, appendVarArgs = F)

setMethodS3("setData", "FactorAnalyses", function(this, X) {
  if (ncol(X) != this$.numCovariates) {
    throw("Number of covariates in model does not match input matrix.")
  }
  this$.X = X
  this$.logLikes = rep(NA, this$getNumModels())
}, appendVarArgs = F)

setMethodS3("getData", "FactorAnalyses", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

setMethodS3("getNumSamples", "FactorAnalyses", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

setMethodS3("parents", "FactorAnalyses", function(this, model) {
  if (model > this$getNumModels() ||
      model <= 0 || length(model) != 1) {
    throw("Invalid input model.")
  }
  if (model == 1) {
    return(numeric(0))
  } else {
    return(model - 1)
  }
}, appendVarArgs = F)

setMethodS3("logLikeMle", "FactorAnalyses", function(this, model, starts = 1) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$getData()
  n = nrow(X)
  m = this$.numCovariates
  S = cov(X)
  k = this$getNumFactorsForModel(model)
  if (k == 0) {
    R = diag(1, m)
    ell.hat <- sum(log(diag(S))) + m
  } else{
    fit <- factanal(X, k, control = list(nstart = starts))
    R <- tcrossprod(fit$loadings) + diag(fit$uniquenesses)
    ell.hat <- sum(log(diag(S))) + log(det(R)) + m
  }
  this$.logLikes[model] = -n / 2 * ell.hat
  return(this$.logLikes[model])
}, appendVarArgs = F)

setMethodS3("learnCoef", "FactorAnalyses", function(this, superModel, subModel) {
  m = this$.numCovariates
  k = this$getNumFactorsForModel(superModel)
  l = this$getNumFactorsForModel(subModel)
  return( list(lambda = 1 / 4 * ((k + 2) * m + l * (m - k + 1)), m = 1))
}, appendVarArgs = F)

setMethodS3("getDimension", "FactorAnalyses", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
