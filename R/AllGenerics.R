### -------------------------------------------------------------------------
### aldex.clr
###

setGeneric("getMonteCarloInstances", function(.object) standardGeneric("getMonteCarloInstances"))

setGeneric("getDirichletInstances", function(.object) standardGeneric("getDirichletInstances"))

setGeneric("getSampleIDs", function(.object) standardGeneric("getSampleIDs"))

setGeneric("numFeatures", function(.object) standardGeneric("numFeatures"))

setGeneric("getFeatures", function(.object) standardGeneric("getFeatures"))

setGeneric("numMCInstances", function(.object) standardGeneric("numMCInstances"))

setGeneric("getFeatureNames", function(.object) standardGeneric("getFeatureNames"))

setGeneric("getReads", function(.object) standardGeneric("getReads"))

setGeneric("numConditions", function(.object) standardGeneric("numConditions"))

setGeneric("getConditions", function(.object) standardGeneric("getConditions"))

setGeneric("getDirichletReplicate", function(.object, i=-1) standardGeneric("getDirichletReplicate"))

setGeneric("getDirichletSample", function(.object, i=-1) standardGeneric("getDirichletSample"))

setGeneric("getMonteCarloReplicate", function(.object, i=-1) standardGeneric("getMonteCarloReplicate"))

setGeneric("getMonteCarloSample", function(.object, i=-1) standardGeneric("getMonteCarloSample"))

setGeneric("getDenom", function(.object) standardGeneric("getDenom"))

setGeneric("getScaleSamples", function(.object) standardGeneric("getScaleSamples"))

setGeneric("getPrior", function(.object) standardGeneric("getPrior"))

setGeneric("aldex.clr", function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, gamma = NULL, prior=0.5) standardGeneric("aldex.clr"), signature=c("reads"))
