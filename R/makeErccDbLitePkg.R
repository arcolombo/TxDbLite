#' create a ErccDbLite package from sqlite file (usually via repDbLiteFromFasta)
#' @param erccdblitefile   the sqlite filename
#' @param author          whose fault this is
#' @param email           email address 
#' @param version         version of the package (default is "1.0")
#' @param destDir         where to put the new package directory (".")
#' @param ...             additional parameters
#' @importFrom Biobase createPackage
#' @importFrom ensembldb metadata
#' @return the name of the package 
#' 
#' @export
makeErccDbLitePkg <- function(erccdblitefile, author="Nobody", email="dev@null.com", version="1.0", destDir=".", ...) {
  
  stopifnot(class(erccdblitefile) == "character")
  erccdb <- ErccDbLite(x=erccdblitefile)
  md <- metadata(erccdb)
  fetchMeta <- function(x) md[x, "value"]
  pkg <- fetchMeta("package_name")

  maintainer <- paste0(author, " <", email, ">")
  template_path <- system.file("erccdblite", package="TxDbLite")
  release_date <- fetchMeta("creation_time")
  symvals <- list(
    PKGTITLE="ERCC annotation package",
    PKGDESCRIPTION="lightweight ERCC spike-in annotations",
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC="Artistic-2.0",
    ORGANISM="N/A", 
    SPECIES="N/A", 
    PROVIDER="ERCC/NIST",
    PROVIDERVERSION="N/A",
    RELEASEDATE=release_date,
    SOURCEURL="http://www.nist.gov/mml/bbd/ercc.cfm",
    ORGANISMBIOCVIEW="ERCC",

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
  db_path <- file.path(destDir, pkg, "inst", "extdata", erccdblitefile)
  file.copy(erccdblitefile, to=db_path)
  return(pkg)
}
