#' map transcripts or genes to Reactome pathways, using a cache if possible
#'
#' @param     IDs      the identifiers to map
#' @param     type     transcript- or gene-level mappings? (default: transcript)
#' @param     species  what species are the IDs from? (default:Homo sapiens)
#' @param     build    what build the IDs are from? (default:84)
#' @param     asNames  boolean, convert the terms to full pathway names? (FALSE)
#' @param     asUrls   boolean, convert terms to Reactome website URLs? (FALSE)
#' @param     useCache use the cache? (TRUE; do not change unless maintaining)
#' 
#' @return    a list of mappings for the supplied IDs
#' 
#' @examples  mapToReactome("ENSG00000234258", type="gene")
#'            mapToReactome("ENST00000387459", asNames=TRUE)
#'
#' @export
mapToReactome <- function(IDs, 
                          type=c("transcript","gene"),
                          species="Homo sapiens", 
                          build=84,
                          asNames=FALSE, 
                          asUrls=FALSE,
                          useCache=TRUE) { 

  if (asNames & asUrls) stop("You can choose AT MOST one of asNames and asUrls")

  abbreviation <- "HSA" # default
  if (species != "Homo sapiens") {
    abbreviations <- TxDbLite:::getSupportedAbbreviations("reactome")
    names(abbreviations) <- tolower(names(abbreviations))
    joined <- tolower(gsub(" ", "_", gsub("\\.", " ", species)))
    if (!joined %in% names(abbreviations)) stop("Unsupported species.")
    else abbreviation <- tolower(abbreviations[joined])
  }
  orgDetails <- getOrgDetails(species)
  
  type <- match.arg(type)  # now screen for this type
  # FIXME: will fail for C elegans at transcript level 
  IDs <- grep(with(orgDetails, switch(type, transcript=txpre, gene=gxpre)),
              IDs, value=TRUE) # only keep IDs that match & can possibly map

  cache <- getReactomeCache(species=species, build=build)
  res <- cache[intersect(IDs, names(cache))] 
  if (asNames == TRUE) {
    if (!exists("reactomePathways")) data(reactomePathways, package="TxDbLite")
    res <- lapply(res, function(x) reactomePathways[x])
  } else if (asUrls == TRUE) {
    res <- lapply(res, getReactomeUrl)
  }
  return(res)
}
