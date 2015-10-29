#' find any duplicate seqnames in FASTA files BEFORE creating an index...
#' FIXME: remove them, write out a de-duped merged FASTA, and index that.
#' 
#' @param  ...        the FASTA file names (may be compressed, doesn't matter)
#'
#' @return data.frame of duplicate seqnames and fasta filenames, else NULL
#'
#' @import Rsamtools
#' @import Biostrings
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
    duped$allIdentical<-NA 
    duped$ID<-rownames(duped)#unique IDs
    dupeSeqs <- DNAStringSet(apply(duped, 1, getDupeSeq))
    duped<-.determineIdentical(duped,dupeSeqs)
    #creates a list each entry has duplicated sequence under the sequence name
    splitDupe<-split(dupeSeqs,duped$seqname)[duped$seqnames]
      .printIdentical(splitDupe)
     return(duped)
  }

  # if no dupes,
  return(NULL)
}



.determineIdentical<-function(duped,dupeSeqs)  {   
    
    for(i in 1:nrow(duped)){
         id<-duped$ID[duped$seqnames %in% duped$seqnames[i]]
         groupID<-dupeSeqs[which(names(dupeSeqs)==id)]
          #need to run toString and check seed match
          #output is identical Y or N
             seed<-groupID[1]
             for ( g in 1:length(groupID)){
                 if(setequal(seed,groupID[g])==TRUE){
                 duped[id,which(colnames(duped)=="allIdentical")][g]<-"Y"   
                  }#{{{ if seed is true
                           
                 if(setequal(seed,groupID[g])==FALSE){
                 duped[id,which(colnames(duped)=="allIdentical")][g]<-"N"
                   }#{{{ if seed FALSE
 
             }#{{{ for group ID
                
              
   } # for
     return(duped)
}#{{{ main
 

.printIdentical<-function(splitDupe){
df<-sapply(splitDupe,function(x) length(x))
  for(j in 1:length(df)){  
   for( i in 1:df[j]){
 stopifnot(lapply(splitDupe[j],function(x) setequal(x[1],x[i]))==TRUE)
 message(paste0("The duplicated repeats ", names(df)[j]  ," at row number ", names(splitDupe[[j]])[i], " have identical sequences: ",lapply(splitDupe[j],function(x) setequal(x[1],x[i]))==TRUE))
  }# for i 
}#for j

}#{{{ printIdentical main
