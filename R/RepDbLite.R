#' Constructor for an RepDbLite object.
#'
#' @param x the name of the sqlite file
#' @param ... additional arguments related to Repetitive annotation database.
#' @return a RepDbLite object.
#' 
#' @export
RepDbLite <- function (x, ...) {
  repdb <- TxDbLite(x, ...)
  class(repdb) <- "RepDbLite"
  return(repdb)
}
