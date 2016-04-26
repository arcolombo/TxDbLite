#' Pulling Reactome mappings through biomaRt takes time. So we try not to do it.
#' More specifically, try to only pull updates, and then only for new builds.
#'
#' @param   species         species for which to get mappings (Homo sapiens)
#' @param   build           numeric ENSEMBL build (84 by default) 
#'
#' @return  a list of REACTOME terms mapped to ENSEMBL IDs (both tx and gene)
#'
#' @import  utils 
#' 
#' @export
getReactomeCache <- function(species="Homo sapiens", build=84) { 

  if (!is.numeric(build)) stop("The build identifier must be numeric.") 
  underscored <- gsub(" ", "_", gsub("\\.", " ", species))
  speciesKeys <- getSupportedAbbreviations("reactome")
  if (!underscored %in% names(speciesKeys)) stop("Unsupported species name.")
  else whichCache <- speciesKeys[underscored] 

  # this step requires some bootstrapping; maintain for Hsapiens, Mmus, Dmel
  mostRecentBuild <- getMostRecentCacheBuild() 
  if (build > mostRecentBuild) {
    message("Cached ENSEMBL build for ", species, " is ", mostRecentBuild, ".")
    message("The build you have requested is ", build, ".")
    message("Please email the maintainer, ", utils::maintainer("TxDbLite"))
    message("and tell them to update the TxDbLite Reactome caches.")
  }

  if (!exists("reactomeCache")) data(reactomeCache, package="TxDbLite")
  return(reactomeCache[[whichCache]])

} 
