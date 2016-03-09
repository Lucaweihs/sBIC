# This set of code provides functions to compute the singular BIC for a set of
# partially ordered models. The partially ordered set of models is represented
# by a directed acyclic graph (DAG). Currently, DAGs are implemented by using the
# the igraph package.

# Any partially order set can be represented as a DAG. In this representation, each node is a model.
# The root node is the minimal element.
# There can be more than one minimal element, and hence there can be more than one root node.
# A model j \leq i iff i can be reached from j

# The singular BIC is computed recursively. This is done by traversing the graph and computing the BIC for each node.
# When traversing the graph, we need to visit the graph in a preorder fashion i.e. visit each node before it's children.
# This can be done by ordering the nodes by a topological ordering of the nodes.

# To make the code as generic as possible, the following requirements must be met:
# The graph that represents the topological ordering of the models must be created by a function
# The computation of mle and the learning coefficients must be handeled independently.

# Currently to implement any new model family, three new functions must be written:
# Let us say the model family is gaussian mixture models and we will represent it by the string "gaussian.mixture"
# The three functions must be:

# Why do we need the model object?
# Each model family will have different parameters and we need to pass these to different functions
# One way to solve this is to create lists of parameters and pass that generic list
# The burden is now on the function that creates the parameter list to make sure it creates
# the right list, so that the calling function can read it.
# for computing the mle, there will be parameters that are fixed for all the models
# for example, in the reduced rank case, the matrix C is to be computed only once, but is used
# by all the models. Such computations will be saved in a list for each model family.

# Computing the sBIC Iteratively

# note that igraph starts labeling the nodes from 0, but in R
# loops start from 1, so you will see a +1 in many places.
# Need to find a better way to solve this issue.
# list of parents for each edge

sBIC = function(X, mod) {
  numModels = mod$getNumModels()
  topOrder = mod$getTopOrder()
  mod$setData(X)
  n <- nrow(X)

  # go through the vertices of g in the topological order and
  # find the list of reachable nodes (the set {j:j<i})
  reach = vector(mode = "list", length = numModels)
  for (i in topOrder) {
    parents = mod$parents(i)
    for (j in parents) {
      if (length(j) != 0) {
        reach[[i]] <- union(reach[[i]], union(reach[[j]], j))
      }
    }
  }
  if (is.null(reach[[topOrder[1]]]) == FALSE) {
    print(
      'Warning: There is no minimal element according to the toplogical order, check the construction of graph'
    )
  }

  L <- rep(0, numModels)
  p <- mod$getPrior()

  # compute p and l
  # How to represent Lij (or LogL)? Lij is a list of vectors.
  # Each list element Lij[[i]] corressponds to the larger model i
  # For each i, Lij[[i]] is a vector of Lij's for the submodels given by the "reach" list
  # The log L'_{ii} are stored in logLii
  logL <-  vector(mode = "list", length = numModels)
  logLii <- rep(0, numModels)
  Lij <-  vector(mode = "list", length = numModels)
  logLij <-  vector(mode = "list", length = numModels)
  loglike <- rep(0, numModels)

  for (i in topOrder) {
    loglike[i] <- mod$logLikeMle(i)
    lf <- mod$learnCoef(i, i)
    logLii[i] <- loglike[i] - lf$lambda * log(n)
    for (j in reach[[i]]) {
      lf <- mod$learnCoef(i, j)
      logL[[i]][j] <-
        loglike[i] - lf$lambda * log(n) + (lf$m - 1) * log(log(n))
    }
  }
  mn = max(c(unlist(logL), logLii), na.rm = TRUE) # note na.rm=TRUE

  # Normalization
  Lii = exp(logLii - mn)
  for (i in topOrder) {
    Lij[[i]] = exp(logL[[i]] - mn)
  }

  # Compute the values of L iteratively for each model using the topological ordering
  L = rep(0, numModels)
  for (i in topOrder) {
    # the first node is the minimal node.
    # There can be more than one minimal nodes!
    if (is.null(reach[[i]]) == TRUE) {
      L[i] <- Lii[i] # Lij[[i]][1]
    } else{
      a <- p[i]
      b <- -Lii[i] * p[i] + sum(L[reach[[i]]] * p[reach[[i]]])
      c <- -sum(Lij[[i]][reach[[i]]] * L[reach[[i]]] * p[reach[[i]]])
      L[i] <- 1 / (2 * a) * (-b + sqrt(b ^ 2 - 4 * a * c))
    }
  }

  L = pmax(L, 0)
  results = list()
  results$logL = log(L) + mn
  results$L = L
  results$loglike = loglike
  results$regBIC = loglike - (mod$getDimension(1:numModels) / 2) * log(n)
  results$modelPoset = mod
  return(results)
}
