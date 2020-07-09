#' Constructor for an EnsDbLite object.
#'
#' @param x the name of the sqlite file
#' @param ... any additional arguments relating to ENSEMBL annotation database.
#' @return an EnsDbLite object.
#' 
#' @export
EnsDbLite <- function (x, ...) {
  ensdb <- TxDbLite(x, ...)
  class(ensdb) <- "EnsDbLite"
  return(ensdb)
}
