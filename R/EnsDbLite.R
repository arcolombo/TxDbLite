#' Constructor for an EnsDbLite object.
#'
#' @param x the name of the sqlite file
#'
#' @return an EnsDbLite object.
#' 
#' @export
#' 
EnsDbLite <- function (x) {
  ensdb <- TxDbLite(x)
  class(ensdb) <- "EnsDbLite"
  return(ensdb)
}
