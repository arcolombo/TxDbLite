#' helper function for findDupes
#' 
#' @param duperow the row
#' 
#' @return  the sequence for that duplicated entry
#' @importFrom Rsamtools scanFaIndex FaFile
#' @importFrom Biostrings getSeq
#' @importFrom GenomeInfoDb seqnames
#' @export
getDupeSeq <- function(duperow) { 
  seqname <- duperow[1]
  fasta <- duperow[2]
  namedGr <- function(gr) { 
    names(gr) <- seqnames(gr)
    return(gr) 
  }
  faFile <- FaFile(fasta)
  gr <- namedGr(scanFaIndex(faFile))[seqname]
  as.character(getSeq(faFile, gr))
}
