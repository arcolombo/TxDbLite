#' Constructor for an MitxDbLite (miTranscriptome) object.
#'
#' @param x the name of the sqlite file
#'
#' @return an MitxDbLite object.
#' 
#' @export
MitxDbLite <- function (x, ...) {
  mitxdb <- TxDbLite(x, ...)
  class(mitxdb) <- "MitxDbLite"
  return(mitxdb)
}
