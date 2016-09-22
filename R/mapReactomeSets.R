#' this creates a reactomePathway set in terms of gene names from a reactomeSet which is a list in terms of reactomeId per ENSEMBL entry.  *This Function is deprecated* 
#' @param reactomePathways  a cache rda object list with all the reactome Ids listing description of pathways
#' @param reactomeSets  a list of reactome Ids per ENSEMBL gene / transcript set
#' @return a list of reactomePathways per associated gene name

mapReactomeSets<-function(reactomePathways, reactomeSets){

indX<- names(reactomePathways) %in% tx.DF$pathway.name

gn.Url<- getReactomeUrl(gs[,1])
gn.DF<-data.frame(gs,gn.Url, stringsAsFactors=FALSE)


data(reactomePathways,package="TxDbLite")

gindX<-names(reactomePathways) %in% gn.DF$pathway.name
for(i in 1:nrow(gn.DF)) {
ginner<-which(gn.DF$pathway.name[i] == names(reactomePathways[gindX]))
gn.DF$Pathway.Description[i]<-reactomePathways[gindX][[ginner]]
}


}#{{{main

