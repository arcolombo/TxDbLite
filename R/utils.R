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
getFastaStub <- function(fastaFile) { 
  sub("\\.fa$", "", sub("\\.fasta$", "", sub(".gz$", "", fastaFile)))
}


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
createAnnotationPackage <- function(fastaFile, author) {
  if (grepl("GRC", fastaFile) && ## ENSEMBL fasta, or an impostor
      (grepl("cdna", fastaFile) || grep("ncrna", fastaFile))) {
    type <- "EnsDbLite"
    ensdblitefile <- ensDbLiteFromFasta(fastaFile)
    pkg <- makeEnsDbLitePkg(ensdblitefile, author)
    message(pkg, " was created in ", getwd(), "... please install it.")
    return(pkg)
  } else if (grepl("RepBase", fastaFile)) {
    type <- "RepDbLite"
    repdblitefile <- repDbLiteFromFasta(fastaFile)
    pkg <- makeRepDbLitePkg(repdblitefile, author)
    message(pkg, " was created in ", getwd(), "... please install it.")
    return(pkg)
  } else if (grepl("ERCC", fastaFile)) { 
    message("ERCC annotations are avalable in TxDbLite.ERCC.fa, skipping...")
    return(NULL)
  } else {
    message("Don't know how to annotate ", fastaFile, ", skipping...")
    return(NULL)
  }
  setwd(oldwd)
}

