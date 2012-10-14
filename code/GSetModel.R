require(rms)
require(MASS)
require(survival)
require(predictiveModeling)
require(illuminaHumanv3.db)
require(limma)
require(mice)

# if(!("mboost" %in% installed.packages()[,1])) {
#   print("Downloading mboost")
#   download.file('http://cran.r-project.org/src/contrib/mboost_2.1-3.tar.gz', destfile='mboost_2.1-3.tar.gz')
#   print("Installing mboost")
#   install.packages("mboost_2.1-3.tar.gz", repos=NULL)
# }
require(mboost)

setRefClass(Class = "PredictiveModel")

#' gbm model using
#' gene set medians to include gene expression data.
#' @author Thomas Kelder
#' @export
GSetModel = setRefClass(
  Class  = "GSetModel",
  contains = "PredictiveModel",
  fields   = c("modelFit", "exprSets", "selectClinical", "modelType", "iterations"),
  methods  = list(
    # Argument exprSets is a list of gene sets (character vectors, each containing
    # entrez gene ids of the genes belonging to the set)
    initialize = function(exprSets, selectClinical = NULL, modelType = "cox", iterations = 1000) {
      .self$exprSets = exprSets
      .self$selectClinical = selectClinical
      .self$modelType = modelType
      .self$iterations = iterations
      return(.self)
    },
    
    customTrain = function(exprData, copyData, clinicalFeaturesData, clinicalSurvData, ...) {
      trainData = prepareData(exprData, copyData, clinicalFeaturesData)
      
      message("Training model...")
      survTime = clinicalSurvData[,1]
      survStatus = clinicalSurvData[,2]
      
      if(modelType == "cox") {
          ## coxph
          fit = coxph(clinicalSurvData ~., data = as.data.frame(trainData))
      }
      else if(modelType == "glm") {
        fit = glmboost(y = Surv(survTime / 365, survStatus), x = as.matrix(trainData), family = CoxPH(), control = boost_control(mstop = iterations))
      }
      else if(modelType == "gbm") {
        gbm.fit(as.matrix(trainData), clinicalSurvData, distribution="coxph", 
                shrinkage=0.001, n.trees=iterations, interaction.depth=4, bag.fraction=0.8, 
                train.fraction=1, verbose=T)
      }
      .self$modelFit = fit
    },
    
    customPredict = function(exprData, copyData, clinicalFeaturesData) {
      data = prepareData(exprData, copyData, clinicalFeaturesData)
      message("Predicting...")
      if(modelType == "cox") predict(.self$modelFit, as.data.frame(data))
      else if(modelType == "glm") {
        predict(.self$modelFit, as.matrix(data))
      }
      else if(modelType == "gbm") {
        predict(.self$modelFit, as.matrix(data), n.trees = iterations)
      }
    },
    
    prepareData = function(exprData, copyData, clinicalFeaturesData) {
      message("Preparing data...")
      clinFormatted = formatClinical(clinicalFeaturesData)
      clinFormatted = as.matrix(complete(mice(clinFormatted)))
      if(!is.null(.self$selectClinical)) clinFormatted = clinFormatted[, selectClinical]

      trainData = clinFormatted
      
      if(!is.null(.self$exprSets)) {
        edata = exprDataToEntrez(exprData)
        medExpr = medianSets(.self$exprSets, edata)
        trainData = cbind(trainData, medExpr)
      }
      print(colnames(trainData))
      message('train data contains ', ncol(trainData), ' variables')
      return(trainData)
    },
    copy = function() {
      r = GSetModel$new(.self$exprSets)
      r$modelFit = .self$modelFit
      r$exprSets = .self$exprSets
      r$selectClinical = .self$selectClinical
      r$modelType = .self$modelType
      r$iterations = .self$iterations
      return(r)
    }
  )
)
