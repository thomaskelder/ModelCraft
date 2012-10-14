###################
## Load networks ##
###################
library(graphite)
library(devtools)
source_url("https://raw.github.com/thomaskelder/myRScripts/master/network-analysis/igraph.functions.R")
source_url("https://raw.github.com/thomaskelder/myRScripts/master/util/utils.R")

graphiteToIGraph = function() {
  appendPathways = function(l, p, name) {
    message(name)
    pe = lapply(p, convertIdentifiers, 'entrez')
   names(pe) = paste0(names(p), " (", name, ")")
    append(l, pe)
  }
  
  pathways = list()
  pathways = appendPathways(pathways, biocarta, "biocarta")
  pathways = appendPathways(pathways, kegg, "kegg")
  pathways = appendPathways(pathways, reactome, "reactome")
  pathways = appendPathways(pathways, nci, "nci")
  
  graph = graph.empty(directed = T)
  for(pn in names(pathways)) {
    p = pathways[[pn]]
    d = graphite::edges(p)
    if(nrow(d) == 0) next
    
    ## Double edges for undirected
    d = rbind(d, d[d$direction == 'undirected', c(2, 1, 3:ncol(d))])
    
    pgraph = graph.data.frame(d)
    m = list(graph)
    m[[pn]] = pgraph
    graph = mergeGraphs(m, firstAsBase=T)
  }
  
  graph
}

netPath = "data/networks"

files = listFiles(netPath, "gml")
graphs = lapply(files, read.graph, "gml")
graphs = lapply(graphs, as.directed)
graphs = lapply(graphs, function(g) {
  V(g)$name = gsub("L:", "", V(g)$identifier)
  g
})
names(graphs) = files

graphs$graphite = graphiteToIGraph()
graph = mergeGraphs(graphs)
graph = simplify(graph, remove.multiple=F, remove.loops=T)
V(graph)$entrez = V(graph)$name

message("Nodes: ", vcount(graph))
message("Edges: ", ecount(graph))

edata = exprDataToEntrez(trainingData$exprData)
expGraph = induced.subgraph(graph, V(graph)[entrez %in% rownames(edata)])

save(expGraph, graph, file = "data/networks/graph.gxna.RData")