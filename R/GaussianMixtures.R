#' Construct a poset of gaussian mixture models.
#'
#' Creates an object representing a collection of gaussian mixture models which
#' have a natural ordering upon them.
#'
#' @name GaussianMixtures
#' @usage GaussianMixtures(maxNumComponents = 1, learnCoefBound = "")
#' @export GaussianMixtures
#'
#' @param maxNumComponents the maximum number of gaussian components to
#'                         consider in a mixture.
#' @param learnCoefBound which bound on learning coefficients to use, can be
#'                       one of
#'                       \itemize{
#'                        \item{"1"}
#'                        \item{"+"}
#'                       }
#' @param restarts the number of random restarts to perform when computing the
#'        MLE.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("GaussianMixtures",
                 function(maxNumComponents = 1, learnCoefBound = "",
                          restarts = 5000) {
                   numModels = maxNumComponents
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (maxNumComponents == 1) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(1, numModels)
                   for (j in topOrder) {
                     dimension[j] = 3 * j - 1
                   }

                   extend(
                     ModelPoset(),
                     "GaussianMixtures",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .learnCoefBound = learnCoefBound,
                     .restarts = restarts
                   )
                 })

#' @rdname   getTopOrder
#' @name     getTopOrder.GaussianMixtures
#' @export   getTopOrder.GaussianMixtures
setMethodS3("getTopOrder", "GaussianMixtures", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.GaussianMixtures
#' @export   getPrior.GaussianMixtures
setMethodS3("getPrior", "GaussianMixtures", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.GaussianMixtures
#' @export   getNumModels.GaussianMixtures
setMethodS3("getNumModels", "GaussianMixtures", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Set data for the gaussian mixture models.
#'
#' Sets the data to be used by the gaussian mixture models when computing MLEs.
#'
#' @name     setData.GaussianMixtures
#' @export   setData.GaussianMixtures
#'
#' @param this the GaussianMixtures object.
#' @param X the data to be set, should be a numeric vector of observations.
NULL
setMethodS3("setData", "GaussianMixtures", function(this, X) {
  this$.X = as.numeric(X)
  this$.logLikes = rep(NA, this$getNumModels())
}, appendVarArgs = F)

#' @rdname   getData
#' @name     getData.GaussianMixtures
#' @export   getData.GaussianMixtures
setMethodS3("getData", "GaussianMixtures", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.GaussianMixtures
#' @export   getNumSamples.GaussianMixtures
setMethodS3("getNumSamples", "GaussianMixtures", function(this) {
  return(length(this$getData()))
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.GaussianMixtures
#' @export   parents.GaussianMixtures
setMethodS3("parents", "GaussianMixtures", function(this, model) {
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

#' A helper function.
#'
#' A helper function.
#'
#' @name logLikeMleHelper
#' @export logLikeMleHelper
NULL
#' @rdname   logLikeMleHelper
#' @name     logLikeMleHelper.GaussianMixtures
#' @export   logLikeMleHelper.GaussianMixtures
setMethodS3("logLikeMleHelper", "GaussianMixtures", function(this, model,
                                                             restarts) {
  X = this$getData()
  N = length(X)
  numModels = this$getNumModels()

  logLike = mclust::Mclust(X, G = model, model = "V")$loglik
  for (i in 1:restarts) {
    my.z = matrix(rexp(N * model), N, model)
    my.z = my.z / rowSums(my.z)
    temp.fit = mclust::me(modelName = "V", data = X, z = my.z)
    if (!is.na(temp.fit$loglik)) {
      if (temp.fit$loglik > logLike) {
        logLike = temp.fit$loglik
      }
    }
  }
  return(logLike)
}, appendVarArgs = F)

#' @rdname   logLikeMle
#' @name     logLikeMle.GaussianMixtures
#' @export   logLikeMle.GaussianMixtures
setMethodS3("logLikeMle", "GaussianMixtures", function(this, model) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  this$.logLikes[model] = this$logLikeMleHelper(model, this$.restarts)
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   learnCoef
#' @name     learnCoef.GaussianMixtures
#' @export   learnCoef.GaussianMixtures
setMethodS3("learnCoef", "GaussianMixtures", function(this, superModel, subModel) {
  if (this$.learnCoefBound == "") {
    throw(paste("Exact learning coefficients for GaussianMixtures are unknown,",
                "use a bound instead."))
  }

  if (this$.learnCoefBound == "1") {
    lambda = 1 / 2 * ((superModel - 1) + 2 * subModel)
  } else if (this$.learnCoefBound == "+") {
    lambda = 1 / 2 * ((2 * superModel - 1) + subModel)
  } else {
    throw(paste("Bound named: '", this$.learnCoefBound, "' is unknown.",
                sep = ""))
  }

  ## set m to one
  learn.coeff = list(lambda = lambda, m = 1)
  return(learn.coeff)
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.GaussianMixtures
#' @export   getDimension.GaussianMixtures
setMethodS3("getDimension", "GaussianMixtures", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
