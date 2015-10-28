#' Constructor for an ErccDbLite object.
#'
#' @param x the name of the sqlite file
#'
#' @return an ErccDbLite object.
#' 
#' @export
ErccDbLite <- function (x, ...) {
  erccdb <- TxDbLite(x, ...)
  class(erccdb) <- "ErccDbLite"
  return(erccdb)
}
