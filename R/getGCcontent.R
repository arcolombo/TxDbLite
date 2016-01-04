#' Get the GC content of every contig in a FASTA. 
#' 
#' @param fastaFile  the FASTA file.
#'
#' @return a GRanges with a gc_content column. 
#' 
#' @import Rsamtools
#' 
#' @export
getGCcontent <- function(fastaFile) { 
  if (!file.exists(paste0(fastaFile, ".fai"))) indexFa(fastaFile)
  gr <- scanFaIndex(fastaFile)
  mcols(gr)$gc_content <- GCcontent(scanFa(fastaFile)) 
  names(gr) <- seqnames(gr)
  return(gr)
}

#' @describeIn getGCcontent
#'
#' helper function
#' 
#' @param seqs  a DNAStringSet, RNAStringSet, or similar
#'
#' @return a vector of GC proportions
#'
#' @import Biostrings
#' 
#' @export
GCcontent <- function(seqs) {
  rowSums(alphabetFrequency(seqs, as.prob=TRUE)[,c("C","G")])
}
