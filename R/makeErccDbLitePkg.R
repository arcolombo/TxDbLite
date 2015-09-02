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
  erccdb <- ErccDbLite(x=erccdblitefile)
  md <- metadata(erccdb)
  fetchMeta <- function(x) md[x, "value"]
  pkg <- fetchMeta("package_name")
  
 
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
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC="Artistic-2.0",
    ORGANISM=organism, 
    SPECIES=organism, 
    PROVIDER="ERCC",
    PROVIDERVERSION=as.character(ercc_version),
    RELEASEDATE=release_date,
    SOURCEURL=source_url, 
    ORGANISMBIOCVIEW=gsub(" ","_", organism), 
    TXDBOBJNAME=pkg
  )

  require(Biobase)
  createPackage(pkgname=pkg,
                destinationDir=destDir,
                originDir=template_path,
                symbolValues=symvals, 
                ...)
  dir.create(paste(c(destDir, pkg, "inst", "extdata"), 
                   collapse=.Platform$file.sep),
             showWarnings=FALSE, recursive=TRUE)
  db_path <- file.path(destDir, pkg, "inst", "extdata", erccdblitefile)
  file.copy(erccdblitefile, to=db_path)
  return(pkg)

}









