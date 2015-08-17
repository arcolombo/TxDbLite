RepDbLite <- function (x) {
  repdb <- TxDbLite(x)
  class(repdb) <- "RepDbLite"
  return(repdb)
}
