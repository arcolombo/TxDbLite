#' map Reactome functional interactions (FIs) 
#'
#' @param   symbols the identifiers to map
#' @param   type    what type of identifier are these? (gene)
#' 
#' @return  a data.frame mapping IDs to HUGO symbols or transcript names 
#' 
#' @import  biomaRt 
#' 
#' @export
mapFI <- function(species) { 

  if (species != "Homo sapiens") stop("Reactome FI only supports Homo sapiens.")

  pathBase <- system.file("extdata", "", package="TxDbLite")
  reactomeFIdb <- paste0(pathBase, "/FIsInGene_031516_with_annotations.txt") 

  # add in Reactome FI support (gene level)
  reactomeFI <- read.table(reactomeFIdb, head=TRUE, sep="\t") # map by ENSG ?!?
  stop("Reactome FI autoannotation is not fully supported yet")

}
