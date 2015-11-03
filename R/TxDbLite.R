#' Create a TxDbLite object (usually called by subclass constructors)
#' 
#' @param x     name of the sqlite file whence the TxDbLite should be created 
#' @param path  where it lives (default: ".")
#'
#' @return   a TxDbLite object 
#'
#' @export
TxDbLite <- function(x, path=".", ...) { 
  options(useFancyQuotes = FALSE)
  if (!grepl(".sqlite$", x)) x <- paste0(x, ".sqlite")
  x <- paste(path, x, sep=.Platform$file.sep)
  lite <- dbDriver("SQLite")
  con <- dbConnect(lite, dbname = x, flags = SQLITE_RO)
  tables <- dbListTables(con)
  names(tables) <- tables
  Tables <- lapply(tables, function(x) 
    colnames(dbGetQuery(con, paste0("select * from ", x, " limit 1"))))
  new("TxDbLite", con = con, tables = Tables)
}
