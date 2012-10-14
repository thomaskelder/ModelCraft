################################################
## GSEA analysis to select relevant gene sets ##
################################################
library(BCC)
library(GSEABase)
library(multicore)
library(limma)

source("code/functions.R")

edata = exprDataToEntrez(trainingData$exprData)

survTime = trainingData$clinicalSurvData[,1]
survStatus = trainingData$clinicalSurvData[,2]

## Gene set scores for expression data
## Load gene sets
gsetsC2 = getGmt("data/c2.all.v3.0.entrez.gmt", geneIdType = EntrezIdentifier())
gsetsC5 = getGmt("data/c5.all.v3.0.entrez.gmt", geneIdType = EntrezIdentifier())

gsets = append(gsetsC2, gsetsC5)
names(gsets) = c(names(gsetsC2), names(gsetsC5))

## Limma on survival status
design = model.matrix(~0 + factor(trainingData$clinicalSurvData[,2], levels = unique(trainingData$clinicalSurvData[,2])))
colnames(design) = paste("surv", unique(trainingData$clinicalSurvData[,2]), sep='_')
fit = lmFit(edata, design)
contrast.matrix = makeContrasts(surv_1 - surv_0, levels = design)
fit.eb = eBayes(contrasts.fit(fit, contrast.matrix))
tscores = fit.eb$t[,1]

## write lists based on top genes
topgenes = topTable(fit.eb, number=nrow(edata))
topgeneSets = as.list(topgenes[,"ID"])
names(topgeneSets) = topgenes[,"ID"]
save(topgeneSets, file = 'data/topgeneSets.RData')

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

enr = performGSEA(gsets)
enrSig = enr[enr < 1E-5]
gsetSig = lapply(gsets[names(enrSig)], geneIds)
names(gsetSig) = names(enrSig)

## Prune gene sets based on overlap
pruneByOverlap = function(pvals, gsets, maxOvl = 0.5) {
  psort = sort(abs(pvals))
  ssort = gsets[names(psort)]
  include = names(psort)
  for(sn in names(psort)) {
    s1 = ssort[[sn]]
    i = which(names(ssort) == sn)
    if(!(sn %in% include)) next
    ovl = sapply(ssort[(i+1):length(ssort)], function(s2) {
      length(intersect(s1, s2)) / min(c(length(s1), length(s2)))
    })
    ovl[is.na(ovl)] = 0
    include = setdiff(include, names(ovl)[ovl >= maxOvl])
  }
  include
}

gsetInclude = pruneByOverlap(enrSig, gsetSig)
gsetTrain = gsetSig[gsetInclude]

save(gsetTrain, file = "data/gSets.RData")