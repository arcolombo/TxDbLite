#' Constructor for an ErccDbLite object.
#'
#' @param x the name of the sqlite file
#' @param ... additional arguments related to ERCC annotation database.
#' @return an ErccDbLite object.
#' 
#' @export
ErccDbLite <- function (x, ...) {
  erccdb <- TxDbLite(x, ...)
  class(erccdb) <- "ErccDbLite"
  return(erccdb)
}
