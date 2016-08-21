#' (try to) figure out what type of annotation package to create
#' @param fastaFile a fasta file of interest 
#' @export
getAnnotationType <- function(fastaFile) {
  if (grepl("(GRC|Rnor|BDGP|WBcel|R64)", fastaFile) && ## common ENSEMBL genomes
      (grepl("cdna", fastaFile) || 
       grepl("ncrna", fastaFile) || 
       grepl("merged", fastaFile))) {
    type <- "EnsDbLite"
  } else if (grepl("RepBase", fastaFile)) {
    type <- "RepDbLite"
  } else if (grepl("ERCC", fastaFile)) { 
    type <- "ErccDbLite"
  } else {
    type <- NULL
  }
  return(type)
}
