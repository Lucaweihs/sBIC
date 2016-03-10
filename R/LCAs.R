setConstructorS3("LCAs",
                 function(maxNumClasses = 1, numVariables = 2,
                          numStatesForVariables = 2, phi = 1) {
                   numModels = maxNumClasses
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (maxNumClasses == 1) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(1, numModels)
                   K = numeric(numModels)
                   for (j in topOrder) {
                     K[j] = j
                     dimension[j] = (j - 1) + numVariables * (numStatesForVariables - 1) * j
                   }

                   extend(
                     ModelPoset(),
                     "LCAs",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .K = K,
                     .maxNumClasses = maxNumClasses,
                     .numVariables = numVariables,
                     .numStatesForVariables = numStatesForVariables,
                     .phi = phi
                   )
                 })

setMethodS3("getTopOrder", "LCAs", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

setMethodS3("getPrior", "LCAs", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

setMethodS3("getNumModels", "LCAs", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

setMethodS3("setData", "LCAs", function(this, X) {
  if (ncol(X) != this$.numVariables) {
    throw("Input data X has the incorrect number of columns.")
  }
  if (!is.data.frame(X)) {
    X = as.data.frame(X)
  }
  this$.X = X
  this$.logLikes = rep(NA, this$getNumModels())
}, appendVarArgs = F)

setMethodS3("getData", "LCAs", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

setMethodS3("parents", "LCAs", function(this, model) {
  if (model > this$getNumModels() ||
      model < 1 || length(model) != 1) {
    throw("Invalid input model.")
  }
  if (model == 1) {
    return(numeric(0))
  } else {
    return(model - 1)
  }
}, appendVarArgs = F)

setMethodS3("logLikeMle", "LCAs", function(this, model) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$getData()
  f = as.formula(paste("cbind(", paste(names(X), collapse =  ","), ")~1"))
  this$.logLikes[model] = poLCA::poLCA(f, X, nclass = model, nrep = 50,
                                       maxiter = 8000, verbose = FALSE)$llik
  return(this$.logLikes[model])
}, appendVarArgs = F)

setMethodS3("learnCoef", "LCAs", function(this, superModel, subModel) {
  S = superModel
  S0 = subModel
  Y = this$.numStatesForVariables
  N = this$.numVariables

  lambda <- 1 / 2 * min(S0 * ((Y - 1) * N + 1) - 1 + (S - S0) * this$.phi,
                        N * (Y - 1) * S + S0 - 1)
  learn.coeff = list(lambda = lambda, m = 1)
  return(learn.coeff)
}, appendVarArgs = F)

setMethodS3("getDimension", "LCAs", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
