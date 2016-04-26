# downloaded from the Reactome website:

ensemblToReactome <- read.csv("Ensembl2Reactome.txt", sep="\t", 
                             header=FALSE, stringsAsFactors=FALSE)[, 1:2]
names(ensemblToReactome) <- c("ID", "term")

reactomeCache <- split(ensemblToReactome[, 1:2], sapply(ensemblToReactome$term, strpop, "-", 2))[ getSupportedAbbreviations("reactome") ]
reactomeCache <- lapply(reactomeCache, function(x) split(x$term, x$ID))

save(reactomeCache, file="TxDbLite/data/reactomeCache.rda", compress="xz")
