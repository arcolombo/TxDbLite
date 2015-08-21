#' @describeIn utils
#' 
#' get the "stub" of a FASTA filename (no .fa, no .fasta, no .gz)
#' 
#' @param fastaFile the filename
#' 
#' @return a string
#' 
#' @export
#'
getFastaStub <- function(fastaFile) { # {{{
  sub("\\.fa$", "", sub("\\.fasta$", "", sub(".gz$", "", fastaFile)))
} # }}}

#' @describeIn utils
#' 
#' figure out what type of annotation package to create
#' 
#' @export
#'
getAnnotationType <- function(fastaFile) {  # {{{
  if (grepl("GRC", fastaFile) && ## ENSEMBL fasta, or an impostor
      (grepl("cdna", fastaFile) || grep("ncrna", fastaFile))) {
    type <- "EnsDbLite"
  } else if (grepl("RepBase", fastaFile)) {
    type <- "RepDbLite"
  } else if (grepl("ERCC", fastaFile)) { 
    type <- "ErccDbLite"
  } else {
    type <- NULL
  }
  return(type)
} # }}}

#' @describeIn utils
#' 
#' create an annotation package for a FASTA file (try to figure out what kind)
#' 
#' @param fastaFile the filename
#' @param author    the author (with email address)
#' 
#' @return the name of the annotation package, or NULL if uncertain how to do it
#' 
#' @export
#'
createAnnotationPackage <- function(fastaFile, author) { # {{{

  type <- getAnnotationType(fastaFile) 
  if (type == "EnsDbLite") { 
    ensdblitefile <- ensDbLiteFromFasta(fastaFile)
    pkg <- makeEnsDbLitePkg(ensdblitefile, author)
  } else if (type == "RepDbLite") {
    repdblitefile <- repDbLiteFromFasta(fastaFile)
    pkg <- makeRepDbLitePkg(repdblitefile, author)
  } else if (grepl("ERCC", fastaFile)) { 
    erccdblitefile <- erccDbLiteFromFasta(fastaFile)
    pkg <- makeErccDbLitePkg(erccdblitefile, author)
  } else {
    message("Don't know how to annotate ", fastaFile, ", skipping...")
    return(NULL)
  }

  message(pkg, " was created in ", getwd(), "... please install it.")
  return(pkg)

} # }}}
