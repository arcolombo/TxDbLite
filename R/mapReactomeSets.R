#' this creates a reactomePathway set in terms of gene names from a reactomeSet which is a list in terms of reactomeId per ENSEMBL entry.   
#' @param reactomePathways  a cache rda object list with all the reactome Ids listing description of pathways
#' @param reactomeSets  a list of reactome Ids per ENSEMBL gene / transcript set
#' @return a list of reactomePathways per associated gene name

mapReactomeSets<-function(reactomePathways, reactomeSets){

indX<- names(reactomePathways) %in% tx.DF$pathway.name



}#{{{main

