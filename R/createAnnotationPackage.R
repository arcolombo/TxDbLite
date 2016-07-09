#' create an annotation package for a FASTA file (try to figure out what kind)
#' 
#' @param fastaFile the filename
#' 
#' @return the name of the annotation package, or NULL if uncertain how to do it
#' 
#' @import Biobase
#' @import DBI
#' @import RSQLite
#'
#' @export
createAnnotationPackage <- function(fastaFile, ...) { # {{{

  type <- getAnnotationType(fastaFile) 
  if (is.null(type)) {
    message("Couldn't find a known type for ", fastaFile, ", skipping...")
    return(NULL)
  } else if (type == "EnsDbLite") { 
    db <- ensDbLiteFromFasta(fastaFile)
    pkg <- makeEnsDbLitePkg(db, ...)
  } else if (type == "RepDbLite") {
    db <- repDbLiteFromFasta(fastaFile)
    pkg <- makeRepDbLitePkg(db, ...)
  } else if (type == "MitxDbLite") {
    db <- mitxDbLiteFromFasta(fastaFile)
    pkg <- makeMitxDbLitePkg(db, ...) 
  } else if (grepl("ERCC", fastaFile)) { 
    db <- erccDbLiteFromFasta(fastaFile)
    pkg <- makeErccDbLitePkg(db, ...) 
  } else {
    message("Don't know how to annotate ", fastaFile, ", skipping...")
    return(NULL)
  }

  # FIXME: so ghetto it hurts 
  cmd <- paste("R CMD INSTALL", pkg)
  system(cmd)
  return(pkg)

} # }}}
