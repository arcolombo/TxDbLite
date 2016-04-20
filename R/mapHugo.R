#' map genes or transcripts to updated HUGO symbols 
#'
#' @param   IDs     the identifiers to map
#' @param   type    what type of identifier are these? (gene)
#' 
#' @return  a data.frame mapping IDs to HUGO symbols or transcript names 
#' 
#' @import  biomaRt 
#' 
#' @export
mapHugo <- function(IDs, type=c("gene","transcript")) { 

  type <- match.arg(type) 
  idtype <- paste("ensembl", type, "id", sep="_") 
  columns <- list(gene=c(ENSG="ensembl_gene_id", 
                         HUGO="hgnc_symbol"),
                  transcript=c(ENST="ensembl_transcript_id",
                               HUGO="hgnc_transcript_name"))

  ensembl <- useEnsembl(biomart="ensembl", 
                        dataset="hsapiens_gene_ensembl") 
  getBM(attributes=columns[[type]], filters=idtype, values=IDs, mart=ensembl)

}
