#' Pulling Reactome mappings through biomaRt takes time. So we try not to do it.
#' More specifically, try to only pull updates, and then only for new builds. this will subset the reactomeCache.rda to the species of selection and gene-/transcrpit appropriate selection which can be easily made into an enrichment gene/transcript set calling makeReactomeEnrichmentSet.R
#'
#' @param   species         species for which to get mappings (Homo sapiens)
#' @param   build           numeric ENSEMBL build (84 by default) 
#' @param   byLevel         character specifying gene or transcript level to return
#' @return  a list of REACTOME terms mapped to ENSEMBL IDs (both tx and gene)
#'
#' @import  utils 
#' 
#' @export
getReactomeCache <- function(species=c("Homo sapiens","Mus musculus"), build=84, byLevel=c("gene","transcript")) { 
  speciesElect<-match.arg(species,c("Homo sapiens", "Mus musculus")) #for right now
  levelElect<-match.arg(byLevel,c("gene","transcript"))
 
  if (!is.numeric(build)) stop("The build identifier must be numeric.") 
  underscored <- gsub(" ", "_", gsub("\\.", " ", species))
  speciesKeys <- getSupportedAbbreviations("reactome")
  if (!underscored %in% names(speciesKeys)) stop("Unsupported species name.")
  else whichCache <- speciesKeys[underscored] 

  # this step requires some bootstrapping; maintain for Hsapiens, Mmus, Dmel
  mostRecentBuild <- getMostRecentCacheBuild() 
  if (build > mostRecentBuild) {
    message("Cached ENSEMBL build for ", species, " is ", mostRecentBuild, ".")
    message("The build you have requested is ", build, ".")
    message("Please email the maintainer, ", utils::maintainer("TxDbLite"))
    message("and tell them to update the TxDbLite Reactome caches.")
  }

  if (!exists("reactomeCache")) data(reactomeCache, package="TxDbLite")
   
  #now subset based on species selection
   specReactome<-reactomeCache[grep(species,reactomeCache$species),]
   if(nrow(specReactome)==0) { 
     message("it does not appear that the reactome cache has your species...")
    }
  
   #subset based on level selection
    if(speciesElect=="Homo sapiens") {
        if(levelElect=="gene") {
     #grep out human ENSG
        specReactome<- specReactome[grep("^ENSG",specReactome$ensemblID),]
         }
        if(levelElect=="transcript") {
      #grep out human ENST
        specReactome<- specReactome[grep("^ENST",specReactome$ensemblID),]
        }
   } #human level select

  
    if(speciesElect=="Mus musculus") {
        if(levelElect=="gene") {
     #grep out human ENSMUSG
      specReactome<- specReactome[grep("^ENSMUSG",specReactome$ensemblID),]
         }
        if(levelElect=="transcript") {
      #grep out human ENSMUST
      specReactome<-specReactome[grep("^ENSMUST",specReactome$ensemblID),]   
        
         }
     } #mouse

#FIX ME: the rda must be creating into a SQLite library

 # return(reactomeCache[[whichCache]])
  return(specReactome)
} 
