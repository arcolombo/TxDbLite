#' map genes or transcripts to updated HUGO symbols 
#'
#' @param   IDs       the identifiers to map
#' @param   type      what type of identifier are these? (transcript)
#' @param   useCache  boolean: use the cache? (TRUE, and should stay that way)
#' 
#' @return  a list mapping IDs to HUGO symbols or HUGO transcript names 
#' 
#' @examples mapHugo("ENST00000003084")
#' 
#' @export
mapHugo <- function(IDs, useCache=TRUE) { 
  cache <- getHugoCache(useCache=useCache)
  cache[intersect(IDs, names(cache))]
}
