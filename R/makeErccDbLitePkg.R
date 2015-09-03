<<<<<<< HEAD
#' create an ErccDbLite packate from an erccDbLite.sqlite (usually via erccDbLiteFromFasta)
#'
#' @param erccdblitefile    the sqlite filename
#' @param author            whose fault this is
#' @param email             adding maintainer's email 
#' @param version           version of the package (default is "1.0")
#' @param destDir           where to put the new package directory (".")

#' @return the name of the package
#'
#' @export
#'
#'

makeErccDbLitePkg<- function(erccdblitefile,author,email,version="1.0",destDir=".", ...) {

stopifnot(class(erccdblitefile) == "character")
=======
#' create a ErccDbLite package from sqlite file (usually via repDbLiteFromFasta)
#'
#' @param erccdblitefile   the sqlite filename
#' @param author          whose fault this is
#' @param email           email address 
#' @param version         version of the package (default is "1.0")
#' @param destDir         where to put the new package directory (".")
#' 
#' @return the name of the package 
#' 
#' @export
#'
makeErccDbLitePkg <- function(erccdblitefile, author="Nobody", email="dev@null.com", version="1.0", destDir=".", ...) {
  
  stopifnot(class(erccdblitefile) == "character")
>>>>>>> 6069cb06980ff47e70507146a92ca8f3e63bb981
  erccdb <- ErccDbLite(x=erccdblitefile)
  md <- metadata(erccdb)
  fetchMeta <- function(x) md[x, "value"]
  pkg <- fetchMeta("package_name")
<<<<<<< HEAD
  
 
   if (grepl("_",pkg)){
   pkg<-gsub("_",".",pkg)
  }
  

if(missing(email)){
   email<-"TommyTrojan@update.com"
  
}
 #type checking the email parameter, and defaulting if it fails
#weak regex for email additions
 if(grepl("<",email)==TRUE || grepl(">",email)==TRUE) {
    email<-gsub("<","",email)
    email<-gsub(">","",email)
 }  

  if(grepl("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}",email)==TRUE){
    message("email has format ASCII@ASCII.ACII ...")
    email<-paste0("<",email)
    email<-paste0(email,">")
  
   } 
   if( grepl("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}", email)==FALSE){
    email<-"<TommyTrojan@update.com>"
   }
  maintainer<-paste(author,email,sep=" ")
  organism <- fetchMeta("organism")
  ercc_version <- md["genome_build", "value"]
  template_path <- system.file("erccdblite", package="TxDbLite")
  source_url <-"http://bio.math.berkeley.edu/kallisto/transcriptomes/"
  release_date <- fetchMeta("creation_time")
  symvals <- list(
    PKGTITLE="Ercc-based annotation package",
    PKGDESCRIPTION="Lightweight Ercc transcript annotations",
=======

  maintainer<-paste(author,email,sep=" ")
  template_path <- system.file("erccdblite", package="TxDbLite")
  release_date <- fetchMeta("creation_time")
  symvals <- list(
    PKGTITLE="ERCC annotation package",
    PKGDESCRIPTION="lightweight ERCC spike-in annotations",
>>>>>>> 6069cb06980ff47e70507146a92ca8f3e63bb981
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC="Artistic-2.0",
<<<<<<< HEAD
    ORGANISM=organism, 
    SPECIES=organism, 
    PROVIDER="ERCC",
    PROVIDERVERSION=as.character(ercc_version),
    RELEASEDATE=release_date,
    SOURCEURL=source_url, 
    ORGANISMBIOCVIEW=gsub(" ","_", organism), 
=======
    ORGANISM="N/A", 
    SPECIES="N/A", 
    PROVIDER="ERCC/NIST",
    PROVIDERVERSION="N/A",
    RELEASEDATE=release_date,
    SOURCEURL="http://www.nist.gov/mml/bbd/ercc.cfm",
    ORGANISMBIOCVIEW="ERCC",
>>>>>>> 6069cb06980ff47e70507146a92ca8f3e63bb981
    TXDBOBJNAME=pkg
  )

  require(Biobase)
  createPackage(pkgname=pkg,
                destinationDir=destDir,
                originDir=template_path,
<<<<<<< HEAD
                symbolValues=symvals, 
=======
                symbolValues=symvals,
>>>>>>> 6069cb06980ff47e70507146a92ca8f3e63bb981
                ...)
  dir.create(paste(c(destDir, pkg, "inst", "extdata"), 
                   collapse=.Platform$file.sep),
             showWarnings=FALSE, recursive=TRUE)
  db_path <- file.path(destDir, pkg, "inst", "extdata", erccdblitefile)
  file.copy(erccdblitefile, to=db_path)
  return(pkg)

}
<<<<<<< HEAD









=======
>>>>>>> 6069cb06980ff47e70507146a92ca8f3e63bb981
