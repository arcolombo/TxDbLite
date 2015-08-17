makeRepDbLitePkg <- function(repdbfile, version, author, destDir=".") { 
  
  stopifnot(class(repdbfile) == "character")
  repdb <- RepDbLite(x=repdbfile)
  con <- dbconn(repdb)
  pkg <- fetchMetadata(con, "package_name")
  organism <- fetchMetadata(con, "Organism")
  ensembl_version <- fetchMetadata(con, "ensembl_version")
  template_path <- system.file("repdblite", package="TxDbLite")
  source_url <- paste0("ftp://ftp.ensembl.org/pub/release-",
                       ensembl_version, "/gtf/", tolower(organism))
  symvals <- list(
    PKGTITLE="Ensembl-based annotation package",
    PKGDESCRIPTION="Lightweight annotation DB derived from Ensembl GTF",
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=author,
    LIC="Artistic-2.0",
    ORGANISM=organism, 
    SPECIES=organism, 
    PROVIDER="Ensembl",
    PROVIDERVERSION=as.character(ensembl_version),
    RELEASEDATE=fetchMetadata(repdb, "Creation time"),
    SOURCEURL=fetchMetadata(repdb, "ensembl_host"),
    ORGANISMBIOCVIEW=gsub(" ","_", organism), 
    TXDBOBJNAME=pkg
  )

  require(Biobase)
  createPackage(pkgname=pkg,
                destinationDir=destDir,
                originDir=template_path,
                symbolValues=symvals)
  dir.create(paste(c(destDir, pkg, "inst", "extdata"), 
                   collapse=.Platform$file.sep),
             showWarnings=FALSE, recursive=TRUE)
  db_path <- file.path(destDir, pkg, "inst", "extdata", repdbfile)
  file.copy(repdbfile, to=db_path)
  return(pkg)

}
