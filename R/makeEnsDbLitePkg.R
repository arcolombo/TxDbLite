#' create an EnsDbLite package from sqlite file (usually via ensDbLiteFromFasta)
#'
#' @param ensdblitefile   the sqlite filename
#' @param author          whose fault this is
#' @param email           adding maintainer's email
#' @param version         version of the package (default is "1.0")
#' @param destDir         where to put the new package directory (".")
#' @param ...             additional parameters
#' @importFrom Biobase createPackage
#' @return the name of the package 
#' 
#' @export
makeEnsDbLitePkg <- function(ensdblitefile, author="Nobody", email="dev@null.com", version="1.0", destDir=".", ...){
  
  stopifnot(class(ensdblitefile) == "character")
  ensdb <- EnsDbLite(x=ensdblitefile)
  md <- metadata(ensdb)
  fetchMeta <- function(x) md[x, "value"]
  pkg <- fetchMeta("package_name")
   maintainer<-paste0(author, "<", email, ">")
  organism <- fetchMeta("organism")
  ensembl_version <- md["genome_build", "value"]
  template_path <- system.file("ensdblite", package="TxDbLite")
  source_url <- paste0("ftp://ftp.ensembl.org/pub/release-", 
                       ensembl_version)
  release_date <- fetchMeta("creation_time")
  symvals <- list(
    PKGTITLE="Ensembl-based annotation package",
    PKGDESCRIPTION="Lightweight Ensembl transcript annotations",
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC="Artistic-2.0",
    ORGANISM=organism, 
    SPECIES=organism, 
    PROVIDER="Ensembl",
    PROVIDERVERSION=as.character(ensembl_version),
    RELEASEDATE=release_date,
    SOURCEURL=source_url, 
    ORGANISMBIOCVIEW=gsub(" ","_", organism), 
    TXDBOBJNAME=pkg
  )


  createPackage(pkgname=pkg,
                destinationDir=destDir,
                originDir=template_path,
                symbolValues=symvals, 
                ...)
  dir.create(paste(c(destDir, pkg, "inst", "extdata"), 
                   collapse=.Platform$file.sep),
             showWarnings=FALSE, recursive=TRUE)
  db_path <- file.path(destDir, pkg, "inst", "extdata", ensdblitefile)
  file.copy(ensdblitefile, to=db_path)
  return(pkg)

}
