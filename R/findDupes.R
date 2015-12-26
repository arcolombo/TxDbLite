#' find any duplicate seqnames in FASTA files BEFORE creating an index...
#' FIXME: remove them, write out a de-duped merged FASTA, and index that.
#' 
#' @param  ...        the FASTA file names (may be compressed, doesn't matter)
#'
#' @return data.frame of duplicate seqnames and fasta filenames, else NULL
#'
#' @import Biostrings
#' @import
#' 
#' @export
findDupes <- function(fastaFiles=NULL) { 
     #repBase names are tab separated.  we defined dupes as matching names and matching sequences
     seqInput<-sapply(fastaFiles,function(x) readDNAStringSet(x)) #read input
      #find duplicated sequences first 
     dupeSeqs<-lapply(seqInput,function(x) x[duplicated(x)])
     tabSplit<-lapply(dupeSeqs,function(x) strsplit(names(x),"\t") )
     filteredSeqs<-lapply(tabSplit,function(x) (sapply(x,"[", c(1))))
     spaceSplit<-lapply(filteredSeqs,function(x) strsplit(x," ") )
     ensemblNames.duplicatedSequences<-lapply(spaceSplit,function(x) (sapply(x,"[",c(1))))
     #ensemblNames with duplicated sequences



     #all names that are duplicated 
     nameSearch<-lapply(seqInput,function(x) strsplit(names(x),"\t") )
     namesFilter<-lapply(nameSearch,function(x) (sapply(x,"[", c(1))) )
     dupeNames<-lapply(namesFilter,function(x) x[duplicated(x)])
     dupeNames<-Filter(length,dupeNames)
     dupeLengths<-sapply(dupeNames,function(x) length(x) ) 
     dupeDF<-as.data.frame(dupeNames)

 
     if(length(dupeLengths)>0) {#if dupe was detected
     #FIX ME: a dupe must match the dupe sequence with dupe names. 
     #find the duplicated names in the list of duplicated sequnces
     #find set of intersection between dupeNames (dupe names) and filteredSeqs(dupe sequences) 
     duplicates<- lapply(ensemblNames.duplicatedSequences,function(x) x[which(x==dupeDF)])
     duplicates<-Filter(length,duplicates)
     duplicateLengths<-lengths(duplicates)
     fileName<-names(duplicates)
     message("there are duplicated sequences: ")
     print(duplicates)
     return(duplicates)
     } #dupeLengths contains a dupe

    if(length(dupeLengths)==0) {
    message("found no duplicated sequence names .... no dupes found")
    return(dupeLengths) #0

    }#no dupes were found 


}
#  else{  #stuff to delete...  chrs does not work with dupes FIX ME

 #   fastaFiles <- as.list(do.call(c, fastaFiles))
  #  allChrs <- do.call(rbind, lapply(fastaFiles, chrs))
 # if (anyDuplicated(allChrs$seqnames)) {
  #  dupes <- allChrs$seqnames[duplicated(allChrs$seqnames)]
  #  duped <- allChrs[which(allChrs$seqnames %in% dupes),]
  #  duped <- duped[order(duped$seqnames),]
  #  duped$allIdentical<-NA 
  #  duped$ID<-rownames(duped)#unique IDs
  #  dupeSeqs <- DNAStringSet(apply(duped, 1, getDupeSeq))
  #  duped<-.determineIdentical(duped,dupeSeqs)
    #creates a list each entry has duplicated sequence under the sequence name
    #splitDupe<-split(dupeSeqs,duped$seqname)[duped$seqnames]
    #  .printIdentical(splitDupe)
   #  return(duped)
  #}#

  # if no dupes,
#  return(NULL)
#}

#} #{{{ main

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

uniqueDupe<-splitDupe[!duplicated(names(splitDupe))]
df<-sapply(uniqueDupe,function(x) length(x))
  for(j in 1:length(df)){  
   for( i in 1:df[j]){
 stopifnot(lapply(uniqueDupe[j],function(x) setequal(x[1],x[i]))==TRUE)
 message(paste0("The duplicated repeats ", names(df)[j]  ," at row number ", names(uniqueDupe[[j]])[i], " have identical sequences: ",lapply(uniqueDupe[j],function(x) setequal(x[1],x[i]))==TRUE))
  }# for i 
}#for j

}#{{{ printIdentical main

