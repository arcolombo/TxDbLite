makeEnsDbLitePkg <- function(ensdbfile, version, author, destDir=".") { 
  
  stopifnot(class(ensdbfile) == "character")
  ensdb <- EnsDbLite(x=ensdbfile)
  con <- dbconn(ensdb)
  pkg <- fetchMetadata(con, "package_name")
  organism <- fetchMetadata(con, "Organism")
  ensembl_version <- fetchMetadata(con, "ensembl_version")
  template_path <- system.file("ensdblite", package="TxDbLite")
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
    RELEASEDATE=fetchMetadata(ensdb, "Creation time"),
    SOURCEURL=fetchMetadata(ensdb, "ensembl_host"),
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
  db_path <- file.path(destDir, pkg, "inst", "extdata", ensdbfile)
  file.copy(ensdbfile, to=db_path)
  return(pkg)

}
