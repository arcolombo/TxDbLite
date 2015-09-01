#' utility functions
#'
utils <- NULL

#' @describeIn utils
#' 
#' get the "stub" of a FASTA filename (no .fa, no .fasta, no .gz)
#' 
#' @export
#'
getFastaStub <- function(fastaFile) { # {{{
  xx <- sub("\\.fa$", "", sub("\\.fasta$", "", sub(".gz$", "", fastaFile)))
  sub("cdna\\.all", "cdna", xx)
} # }}}

#' @describeIn utils
#' 
#' Get the abbreviation for an organism (among those with which we're familiar)
#' 
#' @export
#'
getOrganismAbbreviation <- function(organism) { # {{{

  organism <- sub("\\.", "_", organism)
  abbrs <- c(Danio_rerio="Drerio",
             Homo_sapiens="Hsapiens",                  ## primary
             Mus_musculus="Mmusculus",                 ## primary
             Rattus_norvegicus="Rnorvegicus",          ## primary
             Caenorhabditis_elegans="Celegans",
             Saccharomyces_cerevisiae="Scerevisiae",
             Drosophila_melanogaster="Dmelanogaster")
  return(abbrs[organism])

} # }}}

#' @describeIn utils
#' 
#' get the name of the package/sqlite file for a FASTA-based annotation
#' 
#' @export
#'
getTxDbLiteName <- function(fastaFile) { # {{{

  fastaStub <- getFastaStub(fastaFile)
  type <- getAnnotationType(fastaStub)

  if (type == "ErccDbLite") {
    return("ErccDbLite.ERCC.97") ## autoinstall?
  } else if(!is.null(type)) {
    tokens <- strsplit(fastaStub, "\\.")[[1]]
    organism <- tokens[1] 
    organism <- sub("\\.", "_", ## try & be robust
                    sub("Mmusculus", "Mus_musculus", 
                        sub("Hsapiens", "Homo_sapiens", organism)))
    genomeVersion <- tokens[2]
    if (length(tokens) > 3) {
      version <- tokens[3]
      what <- tokens[4]
    } else { 
      version <- tokens[2]
      what <- tokens[3]
    }
    return(paste(type, organism, version, what, sep="."))
  } else {
    return(NULL)
  }

} # }}}

#' @describeIn utils
#' 
#' figure out what type of annotation package to create
#' 
#' @export
#'
getAnnotationType <- function(fastaFile) {  # {{{
  if (grepl("(GRC|Rnor|BDGP|WBcel|R64)", fastaFile) && ## common ENSEMBL genomes
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
    pkg <- makeEnsDbLitePkg(ensdblitefile, author,email)
  } else if (type == "RepDbLite") {
    repdblitefile <- repDbLiteFromFasta(fastaFile)
    pkg <- makeRepDbLitePkg(repdblitefile, author, email)
  } else if (grepl("ERCC", fastaFile)) { 
    erccdblitefile <- erccDbLiteFromFasta(fastaFile)
    pkg <- makeErccDbLitePkg(erccdblitefile, author, email)
  } else {
    message("Don't know how to annotate ", fastaFile, ", skipping...")
    return(NULL)
  }

  message(pkg, " was created in ", getwd(), "... please install it.")
  return(pkg)

} # }}}
