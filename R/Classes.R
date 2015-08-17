#' A parent class for lightweight Ensembl and RepBase/RepeatMasker annotations
#' 
#' @export
#'
setClass("TxDbLite",
  representation(con="DBIConnection", tables="list"),
  prototype=list(con=NULL, tables=list())
)

setClass("EnsDbLite", contains="TxDbLite")

setClass("RepDbLite", contains="TxDbLite")
