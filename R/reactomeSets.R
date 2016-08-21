#' using the cached Reactome databases, generate gene sets for qusage etc.
#'
#' @param     species  what species are the IDs from? (default:Homo sapiens)
#' @param     type     transcript- or gene-level mappings? (default: transcript)
#' @param     mappedReactome   list of the results output from mapToReactome, downstream to assemble the list into pathways. 
#' @return    gene sets per Reactome term for this species 
#' 
#' @examples  reactomeSets("Homo sapiens") 
#'
#' @export
reactomeSets <- function(species="Homo sapiens", type=c("transcript","gene"), mappedReactome=NULL) {
 
 
   type <- match.arg(type)
   abbreviation <- "HSA" # default
  if (species != "Homo sapiens") {
    abbreviations <- getSupportedAbbreviations("reactome")
    names(abbreviations) <- tolower(names(abbreviations))
    joined <- tolower(gsub(" ", "_", gsub("\\.", " ", species)))
    if (!joined %in% names(abbreviations)) stop("Unsupported species.")
    else abbreviation <- tolower(abbreviations[joined])
  }
  orgDetails <- getOrgDetails(species)
   
  if(is.null(mappedReactome)==FALSE) {
     message(" creating reactome set from mapped reactome list...")
     trueLengths<-nchar(names(mappedReactome))
     mp<-data.frame(term=unlist(mappedReactome))
     mp$ID<-rownames(mp)
     #manually finding the lengths of nchars picked up from dupes
      uniqueValue<-unique(trueLengths)
      idx<-which(nchar(mp$ID)!=uniqueValue)
     if(length(uniqueValue)>1) {
      message("Dave, I'm afraid I found the transcript lengths to be non-uniform, and this causes problems with non-uniform ensembl transcrpit lengths...please input mapToReactome IDs with all same lengths")
      stopifnot(length(uniqueValue)==1)
      }
     if(length(idx) > 0) { 
      message("I detected duplicates...correcting...")
      if(length(uniqueValue)==1) {
      message("I split the reactomeSet, now correcting the ID lengths for duplicates...")
      if( unique(nchar(mp$ID[idx ])) ==(uniqueValue+1)) {
      mp$ID[idx]<-substr((mp$ID[idx]),1,uniqueValue)
      rownames(mp)<-NULL
       }  
     }
   } #if dupes are found when indexing
   reactomeSet<-split(mp$ID,mp$term)
   message("done\n")
   return(reactomeSet[lapply(reactomeSet,length)>0])

  } #if mappedReactome vector is input

  if(is.null(mappedReactome)==TRUE) { 
  build<-getMostRecentCacheBuild()
  data(reactomeCache, package="TxDbLite")
  allTerms=data.frame(term=unlist(reactomeCache[[toupper(abbreviation)]]))
  allTerms$ID<-rownames(allTerms)
  key <- with(orgDetails, switch(type, transcript=txpre, gene=gxpre))
  allTerms <- subset(allTerms, substr(allTerms$ID, 1, nchar(key)) == key)
  

    keyedcacheNames<-subset(names(reactomeCache[[toupper(abbreviation)]]), substr(names(reactomeCache[[toupper(abbreviation)]]), 1, nchar(key)) == key)
   trueLengths<-unique(nchar(keyedcacheNames))
    message("checking if the cache ensembl names have uniform lengths ...")  
   stopifnot(length(trueLengths)==1)
    idx<-which(nchar(allTerms$ID)!=trueLengths)
   

    if(length(idx) > 0 ) { 
     message("I detected duplicates, correcting the dupe names...")
      if(length(trueLengths)==1) { 
      message("I split the reactomeSet, now correcting the ID lengths for duplicates...")
      if(unique(nchar(allTerms$ID[idx]))==(trueLengths+1)) {
       allTerms$ID[idx]<-substr((allTerms$ID[idx]),1,trueLengths)
       }
    }
 }
  #since split was used there will exists dupe names which occur in multiple reactome IDs,  since we know the original lengths of the names, which have to scale back the R defaulted additions of name.dupe1  where "1" or "2" for second dupe ... is added to any dupe.  we can detect these and remove them.  
   reactomeSet<-split(allTerms$ID, allTerms$term)
   message("done\n")
   return(reactomeSet[lapply(reactomeSet,length)>0])

 
} 

 
  

} #{{{
