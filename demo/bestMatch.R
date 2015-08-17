library(stringdist)
bestMatch <- function(query, corpus) corpus[amatch(query, corpus, maxDist=Inf)]
