#' Constructor for an ArrayControlDbLite object.
#'
#' @param x the name of the sqlite file
#' @param ... additional arguments related to ArrayControl annotation database.
#' @return an ErccDbLite object.
#' 
#' @export
ArrayControlDbLite <- function (x, ...) {
  erccdb <- TxDbLite(x, ...)
  class(erccdb) <- "ArrayControlDbLite"
  return(erccdb)
}
