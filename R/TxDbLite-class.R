#' A parent class for lightweight Ensembl and RepBase/RepeatMasker annotations
#' @rdname TxDbLite-class
#' @slot con a DBI Connection
#' @slot tables a list of annotation entries 
#' @export
setClass("TxDbLite",
  representation(con="DBIConnection", tables="list"),
  prototype=list(con=NULL, tables=list())
)
#' A parent class for Ensembl based annotations
#' @rdname TxDbLite-class
#' @export
setClass("EnsDbLite", contains="TxDbLite")

#' Repeat Parent Annotation Class
#' @rdname TxDbLite-class 
#' @export
setClass("RepDbLite", contains="TxDbLite")

#' Ercc Based Annotation Class
#' @rdname TxDbLite-class
#' @export
setClass("ErccDbLite", contains="TxDbLite")
