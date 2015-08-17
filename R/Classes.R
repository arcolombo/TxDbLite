#' A parent class for lightweight Ensembl and RepBase/RepeatMasker annotations
#' 
#' @export
#'
setClass("TxDbLite",
  representation(con="DBIConnection", tables="list"),
  prototype=list(con=NULL, tables=list())
)

#' @export
#'
setClass("EnsDbLite", contains="TxDbLite")

#' @export
#'
setClass("RepDbLite", contains="TxDbLite")
