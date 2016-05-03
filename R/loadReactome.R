#' this function should only ever be used to reload/maintain Reactome caches
#' all the files are as downloaded from the Reactome/Ensembl website but gzipped
#' FIXME: add loading of interaction maps for all of the supported organisms.
#'
#' @param   TxDbLitePath    where the TxDbLite source build root can be found
#'
#' @return  invisibly, a list of the cache and pathway name listings.
#' 
#' @export
loadReactome <- function(TxDbLitePath="~/Dropbox/TxDbLite") { 

  filePath <- system.file("extdata", "reactome", package="TxDbLite")

  ensemblToReactome <- read.csv(paste0(filePath, "/Ensembl2Reactome.txt.gz"), 
                                sep="\t", header=FALSE, 
                                stringsAsFactors=FALSE)[, 1:2]
  names(ensemblToReactome) <- c("ID", "term")
  reactomeOrganisms <- getSupportedAbbreviations("reactome")
  reactomeCache <- split(ensemblToReactome[, 1:2], 
                         sapply(ensemblToReactome$term, 
                                strpop, "-", 2))[reactomeOrganisms]
  reactomeCache <- lapply(reactomeCache, function(x) split(x$term, x$ID))
  save(reactomeCache, 
       file=paste0(TxDbLitePath, "/data/reactomeCache.rda"), 
       compress="xz")

  reactomePathways <- read.csv(paste0(filePath, "/ReactomePathways.txt.gz"),
                               sep="\t", header=FALSE, stringsAsFactors=FALSE)
  reactomePathways <- reactomePathways[!duplicated(reactomePathways[,1]),]
  reactomePathways <- split(reactomePathways[,2], reactomePathways[,1])
  save(reactomePathways, 
       file=paste0(TxDbLitePath, "/data/reactomePathways.rda"),
                   compress="xz")

  res <- list(cache=reactomeCache, 
              pathways=reactomePathways)
  invisible(res) 
}
