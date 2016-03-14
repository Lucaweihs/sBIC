setConstructorS3("GaussianMixtures",
                 function(maxNumComponents = 1, learnCoefBound = "") {
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
                     .learnCoefBound = learnCoefBound
                   )
                 })

setMethodS3("getTopOrder", "GaussianMixtures", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

setMethodS3("getPrior", "GaussianMixtures", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

setMethodS3("getNumModels", "GaussianMixtures", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

setMethodS3("setData", "GaussianMixtures", function(this, X) {
  this$.X = as.numeric(X)
  this$.logLikes = rep(NA, this$getNumModels())
}, appendVarArgs = F)

setMethodS3("getData", "GaussianMixtures", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

setMethodS3("getNumSamples", "GaussianMixtures", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

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

setMethodS3("logLikeMleHelper", "GaussianMixtures", function(this, model,
                                                             restarts) {
  N = length(this$getData())
  numModels = this$getNumModels()

  logLike = mclust::Mclust(this$getData(), G = model, model = "V")$loglik
  for (i in 1:restarts) {
    my.z = matrix(rexp(N * model), N, model)
    my.z = my.z / rowSums(my.z)
    temp.fit = mclust::me(modelName = "V", data = x, z = my.z)
    if (!is.na(temp.fit$loglik)) {
      if (temp.fit$loglik > logLike) {
        logLike = temp.fit$loglik
      }
    }
  }
  return(logLike)
}, appendVarArgs = F)

setMethodS3("logLikeMle", "GaussianMixtures", function(this, model) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  this$.logLikes[model] = this$logLikeMLEHelper(model, 10)
  return(this$.logLikes[model])
}, appendVarArgs = F)

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

setMethodS3("getDimension", "GaussianMixtures", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
