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

## TODO: add a Code object that refers to this script on GitHub