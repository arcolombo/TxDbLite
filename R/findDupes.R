#' find any duplicate seqnames in FASTA files BEFORE creating an index...
#' FIXME: remove them, write out a de-duped merged FASTA, and index that.
#' 
#' @param  ...        the FASTA file names (may be compressed, doesn't matter)
#'
#' @return data.frame of duplicate seqnames and fasta filenames, else NULL
#'
#' @import Rsamtools
#' 
#' @export
findDupes <- function(...) { 

  fastaFiles <- list(...)
  if (is.character(fastaFiles[[1]]) && length(fastaFiles[[1]]) > 1)
    fastaFiles <- as.list(do.call(c, fastaFiles))
   
  allChrs <- do.call(rbind, lapply(fastaFiles, chrs))
  if (anyDuplicated(allChrs$seqnames)) {
    dupes <- allChrs$seqnames[duplicated(allChrs$seqnames)]
    duped <- allChrs[which(allChrs$seqnames %in% dupes),]
    duped <- duped[order(duped$seqnames),] 
    dupeSeqs <- DNAStringSet(apply(duped, 1, getDupeSeq))
    allIdentical <- sapply(split(dupeSeqs, duped$seqname),
                           function(ds) all(ds == ds[1]))
    duped$allIdentical <- allIdentical[duped$seqname] 
    return(duped)
  }

  # if no dupes,
  return(NULL)
}
