#' Create a TxDbLite object (usually called by subclass constructors)
#' 
#' @param x  name of the sqlite file from which the TxDbLite should be created 
#'
#' @return   a TxDbLite object 
#'
#' @export
#'
TxDbLite <- function(x) { 
  options(useFancyQuotes = FALSE)
  lite <- dbDriver("SQLite")
  con <- dbConnect(lite, dbname = x, flags = SQLITE_RO)
  tables <- dbListTables(con)
  names(tables) <- tables
  Tables <- lapply(tables, function(x) 
    colnames(dbGetQuery(con, paste0("select * from ", x, " limit 1"))))
  new("TxDbLite", con = con, tables = Tables)
}
