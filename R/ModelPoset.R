setConstructorS3("ModelPoset",
                 function() { extend(Object(), "Object") },
                 abstract = T)

#' Topological ordering of models.
#'
#' Returns a topological ordering of models in the collection.
#'
#' @name     getTopOrder
#' @export   getTopOrder
#'
#' @param this the model poset object.
NULL
setMethodS3("getTopOrder", "ModelPoset", function(this) {
  throw("getTopOrder is not yet implemented.")
}, appendVarArgs = F)

#' The prior on the models.
#'
#' Returns the unnormalized prior on the collection.
#'
#' @name     getPrior
#' @export   getPrior
#'
#' @param this the model poset object.
NULL
setMethodS3("getPrior", "ModelPoset", function(this) {
  throw("getPrior is not yet implemented.")
}, appendVarArgs = F)

#' Number of models.
#'
#' Returns the number of models in the collection.
#'
#' @name     getNumModels
#' @export   getNumModels
#'
#' @param this the model poset object.
NULL
setMethodS3("getNumModels", "ModelPoset", function(this) {
  throw("getNumModels is not yet implemented.")
}, appendVarArgs = F)

#' Set data for a model poset.
#'
#' Sets the data to be used by a poset of models when computing MLEs.
#'
#' @name     setData
#' @export   setData
#'
#' @param this the model poset object.
#' @param X the data to be set.
NULL
setMethodS3("setData", "ModelPoset", function(this, X) {
  throw("setData is not yet implemented.")
}, appendVarArgs = F)

#' Return the set data.
#'
#' If data has been set for the object using the setData() function
#' then will return that data, otherwise will throw an error.
#'
#' @name     getData
#' @export   getData
#'
#' @param this the object from which to get the data.
NULL
setMethodS3("getData", "ModelPoset", function(this) {
  throw("getData is not yet implemented.")
}, appendVarArgs = F)

#' Number of samples in the set data.
#'
#' If data has been set using the setData method then returns the
#' number of samples in the data. Otherwise throws an error.
#'
#' @name     getNumSamples
#' @export   getNumSamples
#'
#' @param this the BinomialMixtures object.
NULL
setMethodS3("getNumSamples", "ModelPoset", function(this) {
  throw("getData is not yet implemented.")
}, appendVarArgs = F)

#' Parents of a model.
#'
#' Returns the immediate parents of a given model, i.e. those models
#' M that are (in the poset ordering) less than the given model but for
#' which there exists no other model M' such that M < M' < (given model).
#'
#' @name     parents
#' @export   parents
#'
#' @param this the object representing the model poset.
#' @param model the model for which the parents should be found.
NULL
setMethodS3("parents", "ModelPoset", function(this, model) {
  throw("parents is not yet implemented.")
}, appendVarArgs = F)

#' Maximum likelihood for data.
#'
#' Computes the maximum likelihood of a model in the model poset for the
#' data set using the setData command.
#'
#' @name     logLikeMle
#' @export   logLikeMle
#'
#' @param this the object representing the model poset.
#' @param model the model for which the maximum likelihood should be computed.
NULL
setMethodS3("logLikeMle", "ModelPoset", function(this, model) {
  throw("logLike is not yet implemented.")
}, appendVarArgs = F)

#' Learning coefficient
#'
#' Computes the learning coefficient for a model with respect to one of the
#' model's submodels.
#'
#' @name     learnCoef
#' @export   learnCoef
#'
#' @param this the object representing the model poset.
#' @param superModel the larger model of the two input models.
#' @param subModel the submodel of the larger model.
NULL
setMethodS3("learnCoef", "ModelPoset", function(this, superModel, subModel) {
  throw("learnCoef is not yet implemented.")
}, appendVarArgs = F)

#' Model dimension.
#'
#' Computes the dimension of a model in the model poset.
#'
#' @name     getDimension
#' @export   getDimension
#'
#' @param this the object representing the model poset.
#' @param model the model for which the dimension should be computed.
NULL
setMethodS3("getDimension", "ModelPoset", function(this, model) {
  throw("getDim is not yet implemented.")
}, appendVarArgs = F)
