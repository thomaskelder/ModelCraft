require(rms)
require(MASS)
require(survival)
require(predictiveModeling)

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
      ## Make clinical features data numeric
      clinFormatted = apply(clinicalFeaturesData, 2, as.numeric)
      for(cn in colnames(clinFormatted)) {
        if(sum(is.na(clinFormatted[, cn])) == nrow(clinFormatted)) {
          clinFormatted[, cn] = as.numeric(factor(clinicalFeaturesData[, cn]))
        }
      }
      
      ## Annotate expression data to entrez
      edata = exprs(exprData)
      edata.entrez = unlist(mget(rownames(edata), illuminaHumanv3ENTREZID, ifnotfound=NA))
      rownames(edata) = edata.entrez
      edata = edata[!is.na(rownames(edata)),]
      edata = avereps(edata)
      
      ## Calculate set medians for gene expression
      medExpr = lapply(exprSets, function(s) {
        s = intersect(s, rownames(edata))
        apply(edata[s,], 2, median, na.rm = T)
      })
      names(medExpr) = names(exprSets)
      medExpr = as.data.frame(medExpr)
      return(cbind(clinFormatted, medExpr))
    },
    copy = function() {
      r = RsfGSetModel$new(.self$exprSets)
      r$modelFit = .self$modelFit
      return(r)
    }
    )
)
