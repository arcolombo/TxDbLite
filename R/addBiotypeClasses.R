#' use data(repeat_biotypes) and/or data(ensembl_biotypes) to classify biotypes
#'
#' @param txdblite  A txdblite object (should have metadata on its sqlite file)
#' 
#' @export
#'
addBiotypeClasses <- function(txdblite) {

  if (class(txdblite) == "EnsDbLite") {
    data(ensembl_biotypes)
    biotypes <- ensembl_biotypes
  } else if (class(txdblite) == "RepDbLite") {
    data(repeat_biotypes)
    biotypes <- repeat_biotypes
  } else {
    stop(paste("Not sure what to do with an object of class", class(txdblite)))
  }

  md <- metadata(txdblite)
  fetchMeta <- function(x) md[x, "value"]

  stop("This function has not been completed yet")

}
