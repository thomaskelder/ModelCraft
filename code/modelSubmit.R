require(predictiveModeling)
require(BCC)
require(survival)
require(survcomp)
require(MASS)
require(rms)
require(mice)

source('code/functions.R')

# synapseLogin() ### not required if configured for automatic login
trainingData = loadMetabricTrainingData()
load('data/gxnaSets.RData')
load('data/topgeneSets.RData')
load("data/gSets.RData")

submitModel = function(
  model, modelClass, postfix = "", practice = T, sourceFiles = c("code/functions.R"),
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
  submitCompetitionModel(modelName = modelName, trainedModel=m,
                        rFiles=c(sourceFiles, modelClassFile), isPracticeModel=practice)
}

## Variable selection using glmnet
topNum = 5000
gxnaSepSets = separateGeneSets(gxnaSets[1:1000])
formClin = formatClinical(trainingData$clinicalFeaturesData)
formClinImp = as.matrix(complete(mice(formClin)))
formTopGenes = t(edata[names(topgeneSets[1:topNum]),])
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

modelClass = "GSetModel"
source(paste0("code/", modelClass, ".R"))
model = GSetModel

submitModel(model, modelClass, postfix = "gsetSel.", exprSets = selGSet, selectClinical = selClin, modelType = "cox", practice = F)
submitModel(model, modelClass, postfix = "gxnaSel.", exprSets = selGxna, selectClinical = selClin, modelType = "cox", practice = F)
submitModel(model, modelClass, postfix = "topSel.", exprSets = selTopGenes, selectClinical = selClin, modelType = "cox", practice = F)

submitModel(model, modelClass, postfix = "gsetGlm.", exprSets = gsetTrain, modelType = "glm", practice = F)
submitModel(model, modelClass, postfix = "gxnaGlm.", exprSets = gxnaSepSets, modelType = "glm", practice = F)
submitModel(model, modelClass, postfix = "topGlm.", exprSets = topgeneSets[1:topNum], modelType = "glm", practice = F)

submitModel(model, modelClass, postfix = "gsetSelGbm.", exprSets = selGSet, modelType = "gbm", practice = F)
submitModel(model, modelClass, postfix = "gxnaSelGbm.", exprSets = selGxna, modelType = "gbm", practice = F)
submitModel(model, modelClass, postfix = "topSelGbm.", exprSets = selTopGenes, modelType = "gbm", practice = F)

submitModel(model, modelClass, postfix = "gsetGbm.", exprSets = gsetTrain, modelType = "gbm", practice = F)
submitModel(model, modelClass, postfix = "gxnaGbm.", exprSets = gxnaSepSets, modelType = "gbm", practice = F)
submitModel(model, modelClass, postfix = "topGbm.", exprSets = topgeneSets[1:topNum], modelType = "gbm", practice = F)

#submitModel(model, modelClass, postfix = "clinGlm.", exprSets = NULL, modelType = "glm", practice = F)