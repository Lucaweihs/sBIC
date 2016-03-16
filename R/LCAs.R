#' Construct a poset of latent class analysis models.
#'
#' Creates an object representing a collection of LCA models which
#' have a natural ordering upon them.
#'
#' @name LCAs
#' @usage LCAs(maxNumClasses = 1, numVariables = 2,
#'                         numStatesForVariables = 2, phi = 1)
#' @export LCAs
#'
#' @param maxNumClasses the number of classes in the largest LCA model to
#'        considered.
#' @param numVariables the number of observed variables.
#' @param numStatesForVariables the number of states for each observed variable,
#'        at the moment these must all be equal.
#' @param phi parameter determining the prior placed on latent class
#'        probabilities.
#'
#' @return An object representing the collection.
NULL
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

#' @rdname   getTopOrder
#' @name     getTopOrder.LCAs
#' @export   getTopOrder.LCAs
setMethodS3("getTopOrder", "LCAs", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.LCAs
#' @export   getPrior.LCAs
setMethodS3("getPrior", "LCAs", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.LCAs
#' @export   getNumModels.LCAs
setMethodS3("getNumModels", "LCAs", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Set data for the LCA models.
#'
#' Sets the data to be used by the LCA models when computing MLEs.
#'
#' @name     setData.LCAs
#' @export   setData.LCAs
#'
#' @param this the LCAs object.
#' @param X the data to be set, should be an integer valued matrix where each
#'        row represents a single sample from the observed variables.
NULL
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

#' @rdname   getData
#' @name     getData.LCAs
#' @export   getData.LCAs
setMethodS3("getData", "LCAs", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.LCAs
#' @export   getNumSamples.LCAs
setMethodS3("getNumSamples", "LCAs", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.LCAs
#' @export   parents.LCAs
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

#' @rdname   logLikeMle
#' @name     logLikeMle.LCAs
#' @export   logLikeMle.LCAs
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

#' @rdname   learnCoef
#' @name     learnCoef.LCAs
#' @export   learnCoef.LCAs
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

#' @rdname   getDimension
#' @name     getDimension.LCAs
#' @export   getDimension.LCAs
setMethodS3("getDimension", "LCAs", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)

#' Set phi parameter.
#'
#' Set the phi parameter in an LCAs object to a different value.
#'
#' @name     setPhi
#' @export   setPhi
#'
#' @param this the LCAs object.
#' @param phi the new phi value.
NULL
#' @rdname   setPhi
#' @name     setPhi.LCAs
#' @export   setPhi.LCAs
setMethodS3("setPhi", "LCAs", function(this, phi) {
  if (!is.numeric(phi) || length(phi) != 1) {
    throw("Invalid phi value.")
  }
  this$.phi = phi
}, appendVarArgs = F)
