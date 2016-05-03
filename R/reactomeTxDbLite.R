#' this is an exercise making a SQL lite library from reactome text file.
#' @param reactomepathwayTxt path to txt file
#' @param verbose boolean on whether to print out to stdout


reactomeTxDbLite<-function( reactomepathwayTxt,verbose=TRUE) {

stopifnot(file.exists(reactomepathwayTxt)=="TRUE")
if(verbose) cat("Extracting pathway reactome associations...")
txtFile<-read.table(reactomepathwayTxt,sep="\t",stringsAsFactors=FALSE)
colnames(txtFile)<-c("reactomeID","pathWay","species")
if(verbose) cat("done.\n")
if (verbose) cat("Creating the database...") #{{{

CASH<-getReactomeCache(species=species,build=build)

for(i in 1:length(CASH)) { 
for(j in 1:length(CASH[[i]])) {

idx<-which(CASH[[i]][j]==CASH)

  }

}



