#' Using the function getReactomeCache, the species is used to subset the reactomeCache, and the level is used to further subset the reactomeCache into a dataFrame with respect to species and level. FIX ME: we should add this to to the TxDbLite during the last stage of annotation, to prevent over-writing other material however right now I'm just going to use it to create a DF for connecting with enrichment analysis as opposed to connecting to TxDbLite annotations (which would be ideal, but not necessary).  Then the enrichment analysis is done by lifting up the set in terms of reactome IDs using a loop. if this step is robust enough , then the annotation step for ensDbLiteFromFasta can merely call this funciton using the tx or gene id as a key into this data frame.
#' @param reactomeCache path to csv file containing the ensembl, reactome ID and species, reactome All levels csv file with ',' and ';' removed manually, then used awk to create the cache, using extrallAllLevels.sh
#' @param verbose boolean on whether to print out to stdout
#' @param species  character either Human or Mouse for now
#' @param level    character either gene or transcript.
#' @return a data frame with rownames relative to level and size is relative to species selection
#' @export

makeReactomeEnrichmentSet<-function( reactomeCache,verbose=TRUE,species=c("Homo sapiens","Mus musculus",level=c("gene","transcript")) {
#this calls getReactomeCache.R to fetch the species level appropriate reactomeCache.rda for looping into a list ready for enrichmentAnalysis

speciesElect<-match.arg(species,c("Homo sapiens","Mus musculus"))
levelType<-match.arg(level,c("gene","transcript"))


cash<-getReactomeCache(species=speciesElect,build=build,byLevel=levelType )

geneSet<-list()

for(i in 1:nrow(cash)){ 
idx<-which(cash$reactomeID[i]==cash$reactomeID)
geneSet[[i]]<-cash$ensemblID[idx]
names(geneSet)[i]<-cash$reactomeID[i]
 }


return(geneSet)

} #{{{
