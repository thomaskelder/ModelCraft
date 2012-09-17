library(predictiveModeling)
library(BCC)
library(survival)
library(survcomp)
library(MASS)
library(rms)

# synapseLogin() ### not required if configured for automatic login

source("code/gsea.R")
source("code/RsfGSetModel.R")

testModel = function(class, ...) {
  m = class$new(...)
  m$customTrain(
    trainingData$exprData, trainingData$copyData,
    trainingData$clinicalFeaturesData, trainingData$clinicalSurvData
    )
  
  predictions = m$customPredict(
    trainingData$exprData, trainingData$copyData,
    trainingData$clinicalFeaturesData
    )
  trainPerformance = SurvivalModelPerformance$new(predictions, trainingData$clinicalSurvData)
  message("Concordance: ", trainPerformance$getConcordanceIndex()$c.index)
  
  message("Cross-validating model...")
  cvPerformance = crossValidatePredictiveSurvivalModel(
    class$new(...), trainingData$exprData,
    trainingData$copyData,
    trainingData$clinicalFeaturesData,
    trainingData$clinicalSurvData, numFolds = 3
    )
  
  return(list(trainPerformance = trainPerformance, cvPerformance = cvPerformance))
}

testResult = list()
testResult$rsf.c2 = testModel(class = RsfGSetModel, exprSets = selectC2)
testResult$rsf.c5 = testModel(class = RsfGSetModel, exprSets = selectC5)

sapply(testResult, function(x) x$trainPerformance$getConcordanceIndex()$c.index)
sapply(testResult, function(x) x$cvPerformance$trainPerformanceCV$getFoldCIndices())
sapply(testResult, function(x) x$cvPerformance$testPerformanceCV$getFoldCIndices())