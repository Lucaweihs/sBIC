setConstructorS3("BinomialMixtures",
                 function(maxNumComponents = 1, rousseau = F) {
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
                   K = numeric(numModels)
                   for (j in topOrder) {
                     K[j] = j
                     dimension[j] = 2 * j - 1
                   }

                   extend(
                     ModelPoset(),
                     "BinomialMixtures",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .K = K,
                     .rousseau = rousseau
                   )
                 })

setMethodS3("getTopOrder", "BinomialMixtures", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

setMethodS3("getPrior", "BinomialMixtures", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

setMethodS3("getNumModels", "BinomialMixtures", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

setMethodS3("setData", "BinomialMixtures", function(this, X) {
  this$.X = X

  flexmixFit = flexmix::initFlexmix(
    X ~ 1,
    k = this$topOrder,
    model = flexmix::FLXglm(family = "binomial"),
    control = list(minprior = 0),
    nrep = 10,
    verbose = FALSE
  )
  flexmixEll = flexmix::logLik(flexmixFit)

  this$.logLikes = numeric(this$getNumModels())
  for (j in this$getTopOrder()) {
    this$.logLikes[j] = flexmixEll[j]
  }
}, appendVarArgs = F)

setMethodS3("getData", "BinomialMixtures", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

setMethodS3("parents", "BinomialMixtures", function(this, model) {
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

setMethodS3("logLikeMle", "BinomialMixtures", function(this, model) {
  return(this$.logLikes[model])
}, appendVarArgs = F)

setMethodS3("learnCoef", "BinomialMixtures", function(this, superModel, subModel) {
  if (!this$.rousseau) {
    lambda = 1 / 2 * ((superModel - 1) + subModel)
  } else {
    lambda = 1 / 2 * ((2 * subModel - 1) + (superModel - subModel) / 2)
  }

  ## set m to one
  learn.coeff = list(lambda = lambda, m = 1)
  return(learn.coeff)
}, appendVarArgs = F)

setMethodS3("getDimension", "BinomialMixtures", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
