library(predictiveModeling)
library(BCC)
library(survival)
library(survcomp)
library(MASS)
library(rms)
library(mice)

# synapseLogin() ### not required if configured for automatic login

source("code/RsfGSetModel.R")
source("code/GSetModel.R")
load('data/gxnaSets.RData')
load('data/topgeneSets.RData')
load("data/gSets.RData")

testModel = function(class, ...) {
  m <<- class$new(...)
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
  
  return(list(model = m, trainPerformance = trainPerformance, cvPerformance = cvPerformance))
}

printTestResults = function(testResult) {
  message("train performance")
  print(sapply(testResult, function(x) x$trainPerformance$getConcordanceIndex()$c.index))
  message("train CV")
  print(sapply(testResult, function(x) x$cvPerformance$trainPerformanceCV$getFoldCIndices()))
  message("test CV")
  print(sapply(testResult, function(x) x$cvPerformance$testPerformanceCV$getFoldCIndices()))
  print(colMeans(sapply(testResult, function(x) x$cvPerformance$testPerformanceCV$getFoldCIndices())))
}

## Variable selection using glmnet
topNum = 5000
gxnaSepSets = separateGeneSets(gxnaSets[1:1000])
formClin = formatClinical(trainingData$clinicalFeaturesData)
formClinImp = as.matrix(complete(mice(formClin)))
formTopGenes = as.matrix(medianSets(topgeneSets[1:topNum], edata))
formGxnaGenes = t(edata[names(gxnaSepSets),])
formGSets = as.matrix(medianSets(gsetTrain, edata))
varsClin = selectVars(formClinImp)
varsTopGenes = selectVars(formTopGenes)
varsGxna = selectVars(formGxnaGenes)
varsGSet = selectVars(formGSets)
selClin = varsClin$active.index
selTopGenes = topgeneSets[varsTopGenes$active.index]
selGSet = gsetTrain[varsGSet$active.index]
selGxna = gxnaSepSets[varsGxna$active.index]

testResult = list()
testResult$clin = testModel(class = GSetModel, exprSets = NULL, selectClinical = selClin, modelType = "cox")
testResult$clinGlm = testModel(class = GSetModel, exprSets = NULL, modelType = "glm")
testResult$gset = testModel(class = GSetModel, exprSets = gsetTrain, selectClinical = selClin, modelType = "cox")
testResult$gsetSel = testModel(class = GSetModel, exprSets = selGSet, modelType = "cox")
testResult$gxna = testModel(class = GSetModel, exprSets = gxnaSepSets, selectClinical = selClin, modelType = "cox")
testResult$gxnaSel = testModel(class = GSetModel, exprSets = selGxna, selectClinical = selClin, modelType = "cox")
testResult$gsetGlm = testModel(class = GSetModel, exprSets = gsetTrain, modelType = "glm")
testResult$gxnaGlm = testModel(class = GSetModel, exprSets = gxnaSepSets, modelType = "glm")
testResult$topGenesGlm = testModel(class = GSetModel, exprSets = topgeneSets[1:topNum], modelType = "glm")
testResult$topGenes = testModel(class = GSetModel, exprSets = topgeneSets[1:topNum], selectClinical = selClin, modelType = "cox")
testResult$topGenesSel = testModel(class = GSetModel, exprSets = selTopGenes, selectClinical = selClin, modelType = "cox")


#Gene Sets
submitModel(model, modelClass, postfix = "gsetSel.", exprSets = selGSet, selectClinical = selClin, modelType = "cox", practice = F)
submitModel(model, modelClass, postfix = "gsetGlm.", exprSets = gsetTrain, modelType = "glm", practice = F)
#GXNA
submitModel(model, modelClass, postfix = "gxnaSel.", exprSets = selGxna, selectClinical = selClin, modelType = "cox", practice = F)
submitModel(model, modelClass, postfix = "gxnaGlm.", exprSets = gxnaSepSets, modelType = "glm", practice = F)
#Top genes
submitModel(model, modelClass, postfix = "topSel.", exprSets = selTopGenes, selectClinical = selClin, modelType = "cox", practice = F)
submitModel(model, modelClass, postfix = "topGlm.", exprSets = topgeneSets[1:topNum], modelType = "glm", practice = F)


printTestResults(testResult)