#' Function that lets us avoid storing a redundant pile of URLs in a cache
#' 
#' @param term  the Reactome term to look up 
#'
#' @return url  the URL to look it up on the Reactome website
#'
#' @export
getReactomeUrl <- function(term) { 
  return(paste0("http://www.reactome.org/PathwayBrowser/#/", toupper(term)))
}
