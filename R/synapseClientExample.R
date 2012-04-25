#######################################
## Short example on how to load data ##
## from Synapse.                     ##
#######################################
library(synapseClient)
library(Biobase)

## Login to Synapse
synapseLogin()

## Pointer to the ModelCraft project
modelCraftProject = getEntity('syn296536')

## Load the clinical data
clinicalDataEntity = loadEntity('syn274201')
clinicalData = clinicalDataEntity$objects$`tcga[[i]]`

## Load the normalized TCGA BRCA data
tcgaBrcaDataEntity = loadEntity('syn274561')
tcgaBrcaData = tcgaBrcaDataEntity$objects$eset
tcgaBrcaExpr = exprs(tcgaBrcaData)

## Add a Code object that refers to this script on GitHub (a specific revision)
gitScript = Code(list(name = 'Synapse example script on GitHub', parentId = properties(modelCraftProject)$id))
gitScript = createOrGetEntity(gitScript)

gitScript = addGithubTag(gitScript, "https://raw.github.com/thomaskelder/ModelCraft/master/R/synapseClientExample.R")
gitScript = updateEntity(gitScript)