#################################################
## Identify perturbed modules using GXNA       ##
## http://www-stat.stanford.edu/~serban/gxna/  ##
#################################################
library(igraph)
library(devtools)

source_url("https://raw.github.com/thomaskelder/myRScripts/master/network-analysis/igraph.functions.R")
source_url("https://raw.github.com/thomaskelder/myRScripts/master/util/utils.R")

load("data/networks/graph.gxna.RData")

## Functions to deal with GXNA tool
trimNumber = function(x) {
  return(sprintf("%.3f", x)) # keep only three decimals
}

writePheno = function(pheno, fileName) {
  write.table(pheno, file = fileName, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

writeExpression = function(data, probes, fileName) {
  write.table(cbind(probes, apply(data, 2, trimNumber)), file = fileName, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

prepareGXNA = function(data, probes, groups, projectName, basePath = "data/tmp/") {
  expFile = paste0(basePath, paste(projectName, "exp", sep = "."))
  phenoFile = paste0(basePath, paste(projectName, "phe", sep = "."))
  writeExpression(data, probes, expFile)
  writePheno(t(groups), phenoFile)
  list(expFile = expFile, phenoFile = phenoFile)
}

writeProbeMap = function(packageName, arrayName = packageName, basePath = "data/tmp/") {
  library(packageName, character.only=TRUE) # load package
  map = as.list(get(paste(arrayName, "ENTREZID", sep = "")))
  tbl = t(rbind(names(map), map))
  tbl = tbl[!is.na(tbl[,2]),]
  mapFile = paste0(basePath, paste(arrayName, "ann", sep = "."))
  write.table(tbl, file = mapFile,
              row.names=FALSE, col.names=FALSE)
  mapFile
}

tmpPath = 'data/tmp/'
tmpName = 'tmp'
netFile = paste0(tmpPath, 'graph.txt')
write.table(apply(get.edges(expGraph, E(expGraph)), 2, function(x) get.vertex.attribute(expGraph, "name", x)), 
            file = netFile, row.names = F, col.names = F, quote = F)
edata.prb = exprs(trainingData$exprData)
files = prepareGXNA(edata.prb, rownames(edata.prb), trainingData$clinicalSurvData[,2], tmpName)
probeFile = writeProbeMap("illuminaHumanv3.db", "illuminaHumanv3")

## Run GXNA
oldwd = getwd()
setwd(tmpPath)
gxna = paste0(oldwd, '/code/gxna')
pars = ' -nPerms 100 -algoType 1 -flexSize 1 -depth 15 -runDOT 0'
cmd = paste0(gxna, ' -name ', tmpName, ' -mapFile illuminaHumanv3.ann -edgeFile graph.txt ', pars)
system(cmd)

## Read output
lines = readLines(paste0(tmpName, '_000.res'))
gxnaClusters = lapply(lines, function(l) {
  l = strsplit(l, " ")[[1]]
  ## Genes start on 5th column
  ## Last 3 columns are score info
  info = as.numeric(l[(length(l)-2):length(l)])
  ids = l[5:(length(l)-3)]
  list(ids = ids, score = info[1], pvalue = info[3])
})
setwd(oldwd)

gxnaScores = sapply(gxnaClusters, getElement, 'score')
gxnaPvals = sapply(gxnaClusters, getElement, 'pvalue')
hist(gxnaScores)
gxnaSets = lapply(gxnaClusters, getElement, 'ids')
save(gxnaSets, gxnaClusters, file = 'data/gxnaSets.RData')