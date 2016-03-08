setConstructorS3("ModelPoset",
                 function() { extend(Object(), "ModelPoset") },
                 abstract = T)

setMethodS3("getTopOrder", "ModelPoset", function(this) {
  throw("getTopOrder is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("getPrior", "ModelPoset", function(this) {
  throw("getPrior is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("getNumModels", "ModelPoset", function(this) {
  throw("getNumModels is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("setData", "ModelPoset", function(this, X) {
  throw("setData is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("getData", "ModelPoset", function(this) {
  throw("getData is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("parents", "ModelPoset", function(this, model) {
  throw("parents is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("logLike", "ModelPoset", function(this, model) {
  throw("logLike is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("learnCoef", "ModelPoset", function(this, model1, model2) {
  throw("learnCoef is not yet implemented.")
}, appendVarArgs = F)

setMethodS3("getDim", "ModelPoset", function(this, model) {
  throw("getDim is not yet implemented.")
}, appendVarArgs = F)
