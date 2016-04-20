#' helper function for findDupes
#' 
#' @param fastaFile   the name of the FASTA file
#' 
#' return   the index, invisibly
#' 
#' @export
indexIfNoneFound <- function(fastaFile) {
  if (!file.exists(paste0(fastaFile, ".fai"))) {
    message("Indexing ", fastaFile, " to extract sequence names...")
    invisible(indexFa(fastaFile))
  }
}
