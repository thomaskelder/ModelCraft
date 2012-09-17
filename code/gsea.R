################################################
## GSEA analysis to select relevant gene sets ##
################################################
library(BCC)
library(GSEABase)
library(multicore)
library(limma)

# synapseLogin() ### not required if configured for automatic login
trainingData = loadMetabricTrainingData()

edata = trainingData$exprData
edata = exprs(edata)

survTime = trainingData$clinicalSurvData[,1]
survStatus = trainingData$clinicalSurvData[,2]

## Gene set scores for expression data
## Load gene sets
gsetsC2 = getGmt("data/c2.all.v3.0.entrez.gmt", geneIdType = EntrezIdentifier())
gsetsC5 = getGmt("data/c5.all.v3.0.entrez.gmt", geneIdType = EntrezIdentifier())

## Annotate expression data to entrez gene
library(illuminaHumanv3.db)
edata.entrez = unlist(mget(rownames(edata), illuminaHumanv3ENTREZID, ifnotfound=NA))
rownames(edata) = edata.entrez
edata = edata[!is.na(rownames(edata)),]
edata = avereps(edata)

## Limma on survival status
design = model.matrix(~0 + factor(trainingData$clinicalSurvData[,2], levels = unique(trainingData$clinicalSurvData[,2])))
colnames(design) = paste("surv", unique(trainingData$clinicalSurvData[,2]), sep='_')
fit = lmFit(edata, design)
contrast.matrix = makeContrasts(surv_1 - surv_0, levels = design)
fit.eb = eBayes(contrasts.fit(fit, contrast.matrix))
tscores = fit.eb$t[,1]

## geneSetTest
performGSEA = function(gsets) {
  i = 0
  enr = sapply(gsets, function(s) {
    i <<- i + 1
    if(i %% 100 == 0) message(i)
    ids = geneIds(s)
    geneSetTest(names(tscores) %in% ids, tscores, alternative = "either")
  })
  names(enr) = names(gsets)
  enr
}
enrC2 = performGSEA(gsetsC2)
sigC2 = gsetsC2[enrC2 < 1E-10]
selectC2 = lapply(sigC2, geneIds)
names(selectC2) = names(sigC2)

enrC5 = performGSEA(gsetsC5)
sigC5 = gsetsC5[enrC5 < 1E-5]
selectC5 = lapply(sigC5, geneIds)
names(selectC5) = names(sigC5)