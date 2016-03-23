#' Construct a poset of gaussian mixture models.
#'
#' Creates an object representing a collection of gaussian mixture models. There
#' is one model for each fixed number of components from 1 to some specified
#' maximum. In particular each model is identified by a single number
#' specifiying the number of components in the model. Models are naturally
#' ordered by inclusion so that, for example, a model with 2 components comes
#' before a model with 3 or more components.
#'
#' @name GaussianMixtures
#' @usage GaussianMixtures(maxNumComponents = 1, phi = "default",
#'                         restarts = 5000)
#' @export GaussianMixtures
#'
#' @param maxNumComponents the maximum number of gaussian components to
#'                         consider in a mixture.
#' @param phi parameter determining the prior placed on latent class
#'        probabilities.
#' @param restarts the number of random restarts to perform when computing the
#'        MLE.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("GaussianMixtures",
                 function(maxNumComponents = 1, phi = "default",
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

                   if (phi == "default") {
                     phi = (dimension[1] + 1) / 2
                   }

                   extend(
                     MixtureModels(),
                     "GaussianMixtures",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .phi = phi,
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
  this$.mles = rep(list(NA), this$getNumModels())
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

#' @rdname   logLikeMle
#' @name     logLikeMle.GaussianMixtures
#' @export   logLikeMle.GaussianMixtures
setMethodS3("logLikeMle", "GaussianMixtures", function(this, model) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$getData()
  N = length(X)

  fit = mclust::Mclust(X, G = model, model = "V")
  logLike = fit$loglik
  mle = fit$parameters
  for (i in 1:this$.restarts) {
    my.z = matrix(rexp(N * model), N, model)
    my.z = my.z / rowSums(my.z)
    temp.fit = mclust::me(modelName = "V", data = X, z = my.z)
    if (!is.na(temp.fit$loglik)) {
      if (temp.fit$loglik > logLike) {
        logLike = temp.fit$loglik
        mle = temp.fit$parameters
      }
    }
  }
  this$.logLikes[model] = logLike
  this$.mles[[model]] = list(mixWeights = mle$pro, means = mle$mean,
                             vars = mle$vars)
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   mle
#' @name     mle.GaussianMixtures
#' @export   mle.GaussianMixtures
setMethodS3("mle", "GaussianMixtures", function(this, model) {
  if (!is.na(this$.mle[[model]])) {
    return(this$.mle[[model]])
  }
  this$logLikeMle(model)
  return(this$.mle[[model]])
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.GaussianMixtures
#' @export   getDimension.GaussianMixtures
setMethodS3("getDimension", "GaussianMixtures", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
