library(predictiveModeling)
library(BCC)
library(survival)
library(survcomp)
library(MASS)
library(rms)

# synapseLogin() ### not required if configured for automatic login
trainingData = loadMetabricTrainingData()

source("code/gsea.R")

submitModel = function(
  model, modelClass, postfix = "", practice = T, sourceFiles = c(),
  timestamp = format(Sys.time(), "_%Y%m%d_%H.%M.%S"), ...
) {
  
  modelClassFile = paste0("code/", modelClass, ".R")
  source(modelClassFile)
  
  m = model$new(...)
  m$customTrain(
    trainingData$exprData, trainingData$copyData,
    trainingData$clinicalFeaturesData, trainingData$clinicalSurvData
  )
  modelName = paste0(modelClass, "_", postfix, timestamp)
  submitCompetitionModel(modelName = modelName, trainedModel=model,
                        rFiles=c(sourceFiles, modelClassFile), isPracticeModel=practice)
}

modelClass = "RsfGSetModel"
model = RsfGSetModel
submitModel(model, modelClass, postfix = "C5.1e5.", exprSets = selectC5, practice = F)
submitModel(model, modelClass, postfix = "C2.1e10.", exprSets = selectC2, practice = F)
