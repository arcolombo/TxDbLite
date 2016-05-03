#' using the cached Reactome databases, generate gene sets for qusage etc.
#'
#' @param     species  what species are the IDs from? (default:Homo sapiens)
#' @param     type     transcript- or gene-level mappings? (default: transcript)
#' 
#' @return    gene sets per Reactome term for this species 
#' 
#' @examples  reactomeSets("Homo sapiens") 
#'
#' @export
reactomeSets <- function(species="Homo sapiens", type=c("transcript","gene")) {
  
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
  allTerms <- data.frame(term=unlist(reactomeCache[[abbreviation]]))
  allTerms$ID <- rownames(allTerms)
  key <- with(orgDetails, switch(type, transcript=txpre, gene=gxpre))
  allTerms <- subset(allTerms, substr(allTerms$ID, 1, nchar(key)) == key)
  split(allTerms$ID, allTerms$term)

}
