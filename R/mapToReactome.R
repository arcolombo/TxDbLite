#' map genes or transcripts to Reactome pathways 
#'
#' @param   IDs     the identifiers to map
#' @param   type    what type of identifier are these? (transcript)
#' 
#' @return  a data.frame mapping IDs to Reactome IDs and pathways
#' 
#' @examples mapToReactome("ENST00000519619", type="transcript")
#'
#' @import  biomaRt 
#' 
#' @export
mapToReactome <- function(IDs, type=c("transcript", "gene")) { 

  type <- match.arg(type) 
  idtype <- paste("ensembl", type, "id", sep="_") 
  columns <- list(gene=c(ENSG="ensembl_gene_id", 
                         HUGO="hgnc_symbol"),
                  transcript=c(ENST="ensembl_transcript_id",
                               HUGO="hgnc_transcript_name"))

  message("Fetching mappings from biomart...")
  message("(this may take a little while!)")
  ensembl <- useEnsembl(biomart="ensembl", 
                        dataset="hsapiens_gene_ensembl") 
  getBM(attributes=columns[[type]], filters=idtype, values=IDs, mart=ensembl)

}
