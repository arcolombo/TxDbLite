#' using the cached Reactome databases, generate gene sets for qusage etc.
#'
#' @param     species  what species are the IDs from? (default:Homo sapiens)
#' @param     type     transcript- or gene-level mappings? (default: transcript)
#' @param     mappedReactome   list of the results output from mapToReactome, downstream to assemble the list into pathways. 
#' @return    gene sets per Reactome term for this species 
#' 
#' @examples  reactomeSets("Homo sapiens") 
#'
#' @export
reactomeSets <- function(species="Homo sapiens", type=c("transcript","gene"), mappedReactome=NULL) {
 
 
   type <- match.arg(type)
   abbreviation <- "HSA" # default
  if (species != "Homo sapiens") {
    abbreviations <- TxDbLite:::getSupportedAbbreviations("reactome")
    names(abbreviations) <- tolower(names(abbreviations))
    joined <- tolower(gsub(" ", "_", gsub("\\.", " ", species)))
    if (!joined %in% names(abbreviations)) stop("Unsupported species.")
    else abbreviation <- tolower(abbreviations[joined])
  }
  orgDetails <- getOrgDetails(species)
   
  if(is.null(mappedReactome)==FALSE) {
     message(" creating reactome set from mapped reactome list...")
  mp<-data.frame(term=unlist(mappedReactome))
  mp$ID<-rownames(mp)
  reactomeSet<-split(mp$ID,mp$term)
   }

  if(is.null(mappedReactome)==TRUE) { 
  build<-getMostRecentCacheBuild()
  data(reactomeCache, package="TxDbLite")
  allTerms=data.frame(unlist(reactomeCache[[toupper(abbreviation)]]))
  allTerms$ID<-rownames(allTerms)
  key <- with(orgDetails, switch(type, transcript=txpre, gene=gxpre))
  allTerms <- subset(allTerms, substr(allTerms$ID, 1, nchar(key)) == key)
  reactomeSet<-split(allTerms$ID, allTerms$term)
    }
  message("done.\n")
  return(reactomeSet)

}
