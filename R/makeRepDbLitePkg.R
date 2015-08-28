#' create a RepDbLite package from sqlite file (usually via repDbLiteFromFasta)
#'
#' @param repdblitefile   the sqlite filename
#' @param author          whose fault this is
#' @param email           email address 
#' @param version         version of the package (default is "1.0")
#' @param destDir         where to put the new package directory (".")
#' 
#' @return the name of the package 
#' 
#' @export
#'
makeRepDbLitePkg <- function(repdblitefile, author,email, version="1.0", destDir=".", ...) {
  
  stopifnot(class(repdblitefile) == "character")
  repdb <- RepDbLite(x=repdblitefile)
  md <- metadata(repdb)
  fetchMeta <- function(x) md[x, "value"]
  pkg <- fetchMeta("package_name")

pkg <- fetchMeta("package_name")
   if (grepl("_",pkg)){
   pkg<-gsub("_","",pkg)
  }
   if (grepl(".",pkg)){
   pkg<-strsplit(pkg,split='.',fixed=TRUE)[[1]][1]
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
  repbase_version <- as.character(fetchMeta("genome_build"))
  template_path <- system.file("repdblite", package="TxDbLite")
  source_url <- paste0("http://www.girinst.org/server/RepBase/protected/",
                       "RepBase", repbase_version, ".fasta/")
  release_date <- fetchMeta("creation_time")
  symvals <- list(
    PKGTITLE="RepBase-based annotation package",
    PKGDESCRIPTION="lightweight RepBase transcript annotations",
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC="Artistic-2.0",
    ORGANISM=organism, 
    SPECIES=organism, 
    PROVIDER="RepBase",
    PROVIDERVERSION=repbase_version,
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
  db_path <- file.path(destDir, pkg, "inst", "extdata", repdblitefile)
  file.copy(repdblitefile, to=db_path)
  return(pkg)

}
