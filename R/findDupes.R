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
         indexRep<-vector()
    for(i in 1:length(names(dupeSeqs))){
    indexRep[i]<-which(rownames(duped)==names(dupeSeqs))[i]
       }
    #correcting rep seq names in DNAStringSet
    names(dupeSeqs)<-duped$seqnames[indexRep]
    #creates a list each entry has duplicated sequence under the sequence name
    splitDupe<-split(dupeSeqs,duped$seqname)[duped$seqnames]
    allIdentical<-sapply(splitDupe,function(x) all(x==x[1]))
    stopifnot(all(allIdentical=="TRUE")=="TRUE")
     return(duped)
  }

  # if no dupes,
  return(NULL)
}
