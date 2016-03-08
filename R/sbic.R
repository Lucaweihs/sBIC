#This set of code provides functions to compute the singular BIC for a set of
#partially ordered models. The partially ordered set of models is represented
#by a directed acyclic graph (DAG). Currently, DAGs are implemented by using the
#the igraph package.

#Any partially order set can be represented as a DAG. In this representation, each node is a model.
#The root node is the minimal element.
#There can be more than one minimal element, and hence there can be more than one root node.
#A model j \leq i iff i can be reached from j

#The singular BIC is computed recursively. This is done by traversing the graph and computing the BIC for each node.
#When traversing the graph, we need to visit the graph in a preorder fashion i.e. visit each node before it's children.
#This can be done by ordering the nodes by a topological ordering of the nodes.

#To make the code as generic as possible, the following requirements must be met:
#The graph that represents the topological ordering of the models must be created by a function
#The computation of mle and the learning coefficients must be handeled independently.

#Currently to implement any new model family, three new functions must be written:
#Let us say the model family is gaussian mixture models and we will represent it by the string "gaussian.mixture"
#The three functions must be:

#initialize.model.gaussian.mxiture(dat)
#This function creates and intializes a model object that stores all the information about the models
#For example, it contains the relevant data, the partial order of the models, any other parameters related to the model
#dat is a named list supplied to the function, the function must know what the names in the list are!

#compute.log.like.gaussian.mixture(i,model)
#This function must compute the log likelihood of a model i that is indexed in the graph of the model object.
#It must return a list with a named component ell that contains the loglikelihood.

#compute.learn.coeff.gaussian.mixture(i,j,model)
#This function must compute the learning coefficient of a model j when model i is the larger model
#It must return a named list with components lambda and m

#For each model family, these functions are defined in a new file called "model.family.r"

#Why do we need the model object?
#Each model family will have different parameters and we need to pass these to different functions
#One way to solve this is to create lists of parameters and pass that generic list
#The burden is now on the function that creates the parameter list to make sure it creates
#the right list, so that the calling function can read it.
#for computing the mle, there will be parameters that are fixed for all the models
#for example, in the reduced rank case, the matrix C is to be computed only once, but is used
#by all the models. Such computations will be saved in a list for each model family.

#More documentation to be added!


library("igraph", quietly=TRUE);
source('./functions/gaussian.trees.r');

############################################################################
#        INITIALIZE AND  GENERATE THE MODEL OBJECT                         #
############################################################################
initialize.model=function(dat = list(X=X,Y=Y), model.family = "gaussian.trees"){

  a = paste("initialize.model", model.family, sep=".");
  model = do.call(a,list(dat=dat));
  return(model);

}
############################################################################
#                        COMPUTING THE LIKELIHOOD                          #
############################################################################
compute.log.like=function(i,model,model.family){

  a = paste("compute.log.like",model.family,sep=".");
  log.like = do.call(a,list(i=i,model=model));
  return(log.like);

}
############################################################################
#                  COMPUTING THE LEARNING COEFFCIENT                      #
############################################################################

compute.learn.coeff <- function(i,j,model,model.family){

  a = paste("compute.learn.coeff",model.family,sep=".");
  learn.coeff = do.call(a,list(i=i,j=j,model=model));
  return(learn.coeff);

}
############################################################################
#                  COMPUTING THE SINGULAR BIC ITERATIVELY                  #
############################################################################

logAdd = function(...) {
  x = unlist(list(...))
  if(length(x) == 1) {
    return(x)
  }
  x = sort(x)
  y = x[1]
  for(i in 2:length(x)) {
    if(x[i] == -Inf) {
      next
    }
    y = max(y, x[i]) + log1p(exp(-abs(y - x[i])))
  }
  return(y)
}

logAbsSub = function(x,y) {
  if(x < y) {
    z = x
    x = y
    y = z
  }
  z = x + log1p(-exp(y - x))
  return(z)
}

#note that igraph starts labeling the nodes from 0, but in R
#loops start from 1, so you will see a +1 in many places.
#Need to find a better way to solve this issue.
#list of parents for each edge

BICS = function(dat,model.family="gaussian.trees") {
  mod <- initialize.model(dat = dat, model.family = model.family);
  n.models <- vcount(mod$g);
  pa.list <- get.adjlist(mod$g,mode="in");
  n <- dim(dat$X)[1];

  #go through the vertices of g in the topological order and
  #find the list of reachable nodes (the set {j:j<i})
  reach = vector(mode="list",length=n.models);
  for(i in mod$topo.order){
    for(j in pa.list[[i]]){
      if(length(j) != 0){
        reach[[i]] <- union(reach[[i]],union(reach[[j]],j));
      }
    }
  }
  if(is.null(reach[[mod$topo.order[1]]]) == FALSE){
    print('Warning: There is no minimal element according to the toplogical order, check the construction of graph');
  }
  L <- rep(0, n.models)
  if(!is.null(mod$prior) && length(mod$prior) == n.models) {
    p <- mod$prior
  } else {
    p <- rep(1, n.models)
  }

  #compute p and l;
  #How to represent Lij (or LogL)? Lij is a list of vectors.
  #Each list element Lij[[i]] corressponds to the larger model i
  #For each i, Lij[[i]] is a vector of Lij's for the submodels given by the "reach" list
  #The log L'_{ii} are stored in logLii
  logL <-  vector(mode = "list",length = n.models);
  logLii <- rep(0,n.models);
  Lij <-  vector(mode = "list",length = n.models);
  logLij <-  vector(mode = "list",length = n.models);
  loglike <- rep(0,n.models);

  for(i in mod$topo.order){
    loglike[i] <- compute.log.like(i,mod,model.family)$ell;
    lf <- compute.learn.coeff(i,i,mod,model.family);
    logLii[i] <- loglike[i] - lf$lambda*log(n);
    for(j in reach[[i]]){
      lf <- compute.learn.coeff(i,j,mod,model.family);
      logL[[i]][j] <- loglike[i] - lf$lambda*log(n) + (lf$m - 1)*log(log(n))
    }
  }
  mn = max(c(unlist(logL), logLii), na.rm=TRUE); #note na.rm=TRUE

  ###
  # Normalization
  ###
  Lii = exp(logLii - mn);
  for(i in mod$topo.order){
    Lij[[i]] = exp(logL[[i]] - mn);
  }

  #Compute the values of L iteratively for each model using the topological ordering
  L = rep(0, n.models);
  for(i in mod$topo.order){
    #the first node is the minimal node.
    #There can be more than one minimal nodes!
    if(is.null(reach[[i]]) == TRUE){
      L[i] <- Lii[i]; #Lij[[i]][1];
    }else{
      a <- p[i];
      b <- -Lii[i]*p[i] + sum(L[reach[[i]]]*p[reach[[i]]]);
      c <- -sum(Lij[[i]][reach[[i]]]*L[reach[[i]]]*p[reach[[i]]]);
      L[i] <- 1/(2*a)*(-b+sqrt(b^2-4*a*c));
    }
  }

  #   logLii = logLii - mn
  #   for(i in mod$topo.order){
  #      logLij[[i]] = logL[[i]] - mn
  #    }
  #
  #   #Compute the values of L iteratively for each model using the topological ordering
  #   logL = numeric(n.models)
  #   for(i in mod$topo.order){
  #     #the first node is the minimal node.
  #     #There can be more than one minimal nodes!
  #     if(is.null(reach[[i]]) == TRUE){
  #       logL[i] <- logLii[i];
  #     }else{
  #       #logLii[i] = max(c(logLii[i], logLij[[i]][reach[[i]]])
  #       loga <- log(p[i])
  #       logAbsb <- logAbsSub(logAdd(logL[reach[[i]]] + log(p[reach[[i]]])), logLii[i]+loga)
  #       bSign = sign((logLii[i]+loga <= logAdd(logL[reach[[i]]] + log(p[reach[[i]]]))) - 1/2)
  #       logNegc <- logAdd(logLij[[i]][reach[[i]]] + logL[reach[[i]]] + log(p[reach[[i]]]))
  #       #b <- -(exp(logLii[i]+log(a))) + exp(logAdd(logL[reach[[i]]] + log(p[reach[[i]]])));
  #       #c <- -exp(logAdd(logLij[[i]][reach[[i]]] + logL[reach[[i]]] + log(p[reach[[i]]])));
  #       #logL[i] <- -log(2*a) + log(-b+sqrt(exp(logAdd(2*log(abs(b)), log(4) + log(a) + log(-c)))));
  #       if(logAbsb > 1/2*logAdd(2*logAbsb, log(4) + loga + logNegc)) {
  #         logL[i] <- -log(2) - loga + logAbsb + log1p(exp(1/2*logAdd(2*logAbsb, log(4) + loga + logNegc) - logAbsb))
  #       } else {
  #         logL[i] <- -log(2) - loga + (1/2*logAdd(2*logAbsb, log(4) + loga + logNegc)) + log1p(-bSign*exp(logAbsb - (1/2*logAdd(2*logAbsb, log(4) + loga + logNegc))))
  #       }
  #     }
  #   }

  #mod$logL = logL + mn
  #mod$L <- exp(logL + mn)
  L = pmax(L, 0)
  mod$logL = log(L) + mn
  mod$L <- L
  mod$loglike <- loglike
  mod$regBIC <- loglike - (mod$dimension/2)*log(n)
  return(mod)
}


###############################################################
#           COMPUTING THE SINGULAR BIC APPROXIMATION          #
###############################################################

BICSApprox = function(dat,model.family="gaussian.trees") {

  mod <- initialize.model(dat = dat, model.family = model.family);
  n.models <- vcount(mod$g);
  pa.list <- get.adjlist(mod$g,mode="in");
  n <- dim(dat$X)[1];

  #go through the vertices of g in the topological order and
  #find the list of reachable nodes (the set {j:j<i, not exist k, j<k<i})
  reach = pa.list
  if(length(reach[[mod$topo.order[1]]]) != 0){
    print('Warning: There is no minimal element according to the toplogical order, check the construction of graph');
  }
  L <- rep(0,n.models);
  if(!is.null(mod$prior) && length(mod$prior) == n.models) {
    p <- mod$prior
  } else {
    p <- rep(1, n.models)
  }

  #compute p and l;
  #How to represent Lij (or LogL)? Lij is a list of vectors.
  #Each list element Lij[[i]] corressponds to the larger model i
  #For each i, Lij[[i]] is a vector of Lij's for the submodels given by the "reach" list
  #The log L'_{ii} are stored in logLii
  logL <-  vector(mode = "list",length = n.models);
  logLii <- rep(0,n.models);
  Lij <-  vector(mode = "list",length = n.models);
  loglike <- rep(0,n.models);


  for(i in mod$topo.order){
    loglike[i] <- compute.log.like(i,mod,model.family)$ell;
    lf <- compute.learn.coeff(i,i,mod,model.family);
    logLii[i] <- loglike[i] - lf$lambda*log(n);
    for(j in reach[[i]]){
      lf <- compute.learn.coeff(i,j,mod,model.family);
      logL[[i]][j] <- loglike[i] - lf$lambda*log(n) + (lf$m - 1)*log(log(n))
    }
  }
  mn=max(c(unlist(logL),logLii),na.rm=TRUE); #note na.rm=TRUE

  ###
  # Normalization
  ###
  Lii = exp(logLii - mn);
  for(i in mod$topo.order){
    Lij[[i]] = exp(logL[[i]]-mn);
  }

  #Compute the values of L iteratively for each model using the topological ordering
  L = rep(0,n.models);
  for(i in mod$topo.order){
    #the first node is the minimal node.
    #There can be more than one minimal nodes!
    if(length(reach[[i]]) == 0){
      L[i] <- Lii[i]; #Lij[[i]][1];
    }else{
      a <- p[i];
      b <- -Lii[i]*p[i] + sum(Lii[reach[[i]]]*p[reach[[i]]]);
      c <- -sum(Lij[[i]][reach[[i]]]*Lii[reach[[i]]]*p[reach[[i]]]);
      L[i] <-  1/(2*a)*(-b+sqrt(b^2-4*a*c));
    }
  }

  L = pmax(L, 0)
  mod$logL = log(L) + mn
  mod$L <- L
  mod$loglike <- loglike
  mod$regBIC <- loglike - (mod$dimension/2)*log(n);
  return(mod);
}
