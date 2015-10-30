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
  ## repbase FASTAs can have dupes
  sub("\\.humsub", ".merged", xx)
  sub("\\.humrep", ".merged", xx)
} # }}}

#' @describeIn utils
#' 
#' get a column of a data.frame _complete_with_rownames_ as a vector 
#' 
#' @export
#'
colAsVector <- function(x, y) t(x[, y, drop=FALSE])[1,]

#' @describeIn utils 
#' 
#' supported organism abbreviations
#'
#' @export
#' 
getSupportedAbbreviations <- function() { # {{{
  data(supportedOrganismsForTxDbLite, package="TxDbLite")
  colAsVector(supportedOrganismsForTxDbLite, "abbreviation")
} # }}}

#' @describeIn utils
#' 
#' Get the abbreviation for an organism (among those with which we're familiar)
#' 
#' @export
#'
getOrganismAbbreviation <- function(organism) { # {{{
  org <- sub(" ", "_", sub("\\.", "_", organism))
  abbr <- getSupportedAbbreviations()
  if (org %in% abbr) {
    return(abbr)
  } else if (org %in% names(abbr)) {
    return(abbr[org])
  } else { 
    message(organism, " was not found in data(supportedOrganismsForTxDbLite)")
    message("Supported as 'Genus.species', 'Genus_species' or 'Genus species':")
    for (i in names(abbr)) message("  ", i, " (", abbr[i], " in package names)")
    stop(paste(organism, "was not matched... but pull requests are accepted!"))
  }
} # }}}

#' @describeIn utils
#' 
#' helper fn that handles a number of annoying tasks using saved data
#' 
#' @export
#'
getOrgDetails <- function(organism) { # {{{
  abbr <- getSupportedAbbreviations()
  if (organism %in% abbr) organism <- names(abbr)[which(abbr == organism)] 
  data(supportedOrganismsForTxDbLite, package="TxDbLite")
  return(supportedOrganismsForTxDbLite[organism,])
} # }}}

#' @describeIn utils
#' 
#' get the name of the package/sqlite file for a FASTA-based annotation
#' NOTE: as of 1.9.25, this is based on DbType.Org.Version
#'
#' @export
#'
getTxDbLiteName <- function(fastaFile) { # {{{

  fastaStub <- getFastaStub(fastaFile)
  type <- getAnnotationType(fastaStub)

  if (is.null(type)) {
    return(NULL)
  } else if (type == "ErccDbLite") {
    return("ErccDbLite.ERCC.97") ## autoinstall?
  } else if(!is.null(type)) {
    tokens <- strsplit(fastaStub, "\\.")[[1]]
    organism <- tokens[1] 
    organism <- sub("\\.", "_", ## try & be robust
                    sub("Mmusculus", "Mus_musculus", 
                        sub("Hsapiens", "Homo_sapiens", organism)))
    organism <- getOrganismAbbreviation(organism)
    genomeVersion <- tokens[2]
    if (length(tokens) > 3) {
      version <- tokens[3]
      what <- tokens[4]
    } else { 
      version <- tokens[2]
      what <- tokens[3]
    }
    return(gsub("_", "", paste(type, organism, version, sep=".")))
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
      (grepl("cdna", fastaFile) || 
       grepl("ncrna", fastaFile) || 
       grepl("merged", fastaFile))) {
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
#' 
#' @return the name of the annotation package, or NULL if uncertain how to do it
#' 
#' @export
#'
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
  } else if (grepl("ERCC", fastaFile)) { 
    db <- erccDbLiteFromFasta(fastaFile)
    pkg <- makeErccDbLitePkg(db, ...) 
  } else {
    message("Don't know how to annotate ", fastaFile, ", skipping...")
    return(NULL)
  }

  cmd <- paste("R CMD INSTALL", pkg)
  system(cmd)
  return(pkg)

} # }}}
