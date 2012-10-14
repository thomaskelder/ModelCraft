require(glmnet)

selectVars = function(trainData, nit = 1E5) {
  survTime = trainingData$clinicalSurvData[,1]
  survStatus = trainingData$clinicalSurvData[,2]
  cv.fit = cv.glmnet(trainData, Surv(survTime, survStatus), family = "cox", maxit = nit)
  fit = glmnet(trainData, Surv(survTime, survStatus), family = "cox", maxit = nit)
  Coefficients <- coef(fit, s = cv.fit$lambda.min)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  list(active.coef = Active.Coefficients, active.index = Active.Index)
}

separateGeneSets = function(sets) {
  lsets = as.list(unique(unlist(sets)))
  names(lsets) = lsets
  lsets
}

exprDataToEntrez = function(exprData) {
  library(illuminaHumanv3.db)
  ## Annotate expression data to entrez
  edata = exprs(exprData)
  edata.entrez = unlist(mget(rownames(edata), illuminaHumanv3ENTREZID, ifnotfound=NA))
  rownames(edata) = edata.entrez
  edata = edata[!is.na(rownames(edata)),]
  avereps(edata)
}

medianSets = function(exprSets, edata) {
  ## Calculate set medians for gene expression
  medExpr = lapply(exprSets, function(s) {
    s = intersect(s, rownames(edata))
    apply(edata[s,,drop=F], 2, median, na.rm = T)
  })
  names(medExpr) = names(exprSets)
  as.data.frame(medExpr)
}

formatClinical = function(clinicalFeaturesData) {
  clinFormatted = apply(clinicalFeaturesData, 2, as.numeric)
  for(cn in colnames(clinFormatted)) {
    if(sum(is.na(clinFormatted[, cn])) == nrow(clinFormatted)) {
      clinFormatted[, cn] = as.numeric(factor(clinicalFeaturesData[, cn]))
    }
  }
  clinFormatted
}