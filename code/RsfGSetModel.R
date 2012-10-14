require(rms)
require(MASS)
require(survival)
require(predictiveModeling)
require(illuminaHumanv3.db)
require(limma)
source("code/functions.R")

setRefClass(Class = "PredictiveModel")

#' Random Survival Forest model using
#' gene set medians to include gene expression data.
#' @author Thomas Kelder
#' @export
RsfGSetModel = setRefClass(
  Class  = "RsfGSetModel",
  contains = "PredictiveModel",
  fields   = c("modelFit", "exprSets"),
  methods  = list(
    # Argument exprSets is a list of gene sets (character vectors, each containing
    # entrez gene ids of the genes belonging to the set)
    initialize = function(exprSets) {
      .self$exprSets = exprSets
      return(.self)
    },
    
    customTrain = function(exprData, copyData, clinicalFeaturesData, clinicalSurvData, ...) {
      trainData = prepareData(exprData, copyData, clinicalFeaturesData)
      
      message("Training RSF model...")
      survTime = clinicalSurvData[,1]
      survStatus = clinicalSurvData[,2]
      fit = rsf(
        Surv(survTime, survStatus) ~ ., 
        cbind(trainData, survTime, survStatus), 
        na.action = "na.impute", 
        ntree=1000
        ) 
      
      .self$modelFit = fit
    },
    
    customPredict = function(exprData, copyData, clinicalFeaturesData) {
      data = prepareData(exprData, copyData, clinicalFeaturesData)
      message("Predicting...")
      predict(.self$modelFit, data)$mortality
    },
    
    prepareData = function(exprData, copyData, clinicalFeaturesData) {
      message("Preparing data...")
      clinFormatted = formatClinical(clinicalFeaturesData)
      
      edata = exprDataToEntrez(exprData)
      medExpr = medianSets(exprSets, edata)
      return(cbind(clinFormatted, medExpr))
    },
    copy = function() {
      r = RsfGSetModel$new(.self$exprSets)
      r$modelFit = .self$modelFit
      return(r)
    }
    )
)
