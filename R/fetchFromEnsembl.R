#' a wrapper for biomaRt::getBM 
#' 
#' @param cols    attributes to fetch 
#' @param idtype  ID type filter (or a list of filters, more rarely)
#' @param IDs     the IDs to try and fetch (NULL if filters is a list)
#' @param dataset the dataset to use 
#'
#' @return        a data.frame
#' 
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt getBM
#'
#' @export
fetchFromEnsembl <- function(cols, idtype, IDs=NULL, dataset) {
  # this may need changing now and then ... biomart is flakey like that
  ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset=dataset)
  if (is.null(IDs)) {
    getBM(attributes=cols, filters=idtype, mart=ensembl)
  } else {
    getBM(attributes=cols, filters=idtype, values=IDs, mart=ensembl)
  }
}
