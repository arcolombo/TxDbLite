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
  sub("\\.merged", "", xx)
} # }}}

#' @describeIn utils 
#' 
#' supported organism abbreviations
#'
#' @export
#' 
getSupportedAbbreviations <- function() { # {{{
  abbrs <- c(Danio_rerio="Drerio",
             Homo_sapiens="Hsapiens",                  ## primary
             Mus_musculus="Mmusculus",                 ## primary
             Rattus_norvegicus="Rnorvegicus",          ## primary
             Caenorhabditis_elegans="Celegans",
             Saccharomyces_cerevisiae="Scerevisiae",
             Drosophila_melanogaster="Dmelanogaster")
  return(abbrs)
} # }}}

#' @describeIn utils
#' 
#' Get the abbreviation for an organism (among those with which we're familiar)
#' 
#' @export
#'
getOrganismAbbreviation <- function(organism) { # {{{
  organism <- sub("\\.", "_", organism)
  abbr <- getSupportedAbbreviations()
  if (organism %in% abbr) {
    return(abbr)
  } else if (organism %in% names(abbr)) {
    return(abbr[organism])
  } else { 
    stop(paste(organism, "is not recognized... pull requests accepted!")) 
  }
} # }}}

#' @describeIn utils
#' 
#' helper fn that handles a number of annoying tasks
#' 
#' @export
#'
getOrgDetails <- function(organism) { # {{{
  abbr <- getSupportedAbbreviations()
  if (organism %in% names(abbr)) organism <- abbr[organism] 
  if (organism == "Hsapiens") { 
    package <- "org.Hs.eg.db"
    keytype <- "ENSEMBL"
    symbol <- "SYMBOL"
    txpre <- "ENST"
    gxpre <- "ENSG"
  } else if (organism == "Drerio") {
    package <- "org.Dr.eg.db"
    keytype <- "ENSEMBL"
    symbol <- "SYMBOL"
    txpre <- "ENSDART"
    gxpre <- "ENSDARG"
  } else if (organism == "Mmusculus") {
    package <- "org.Mm.eg.db"
    keytype <- "ENSEMBL"
    symbol <- "SYMBOL"
    txpre <- "ENSMUST"
    gxpre <- "ENSMUSG"
  } else if (organism == "Rnorvegicus") {
    package <- "org.Rn.eg.db"
    keytype <- "ENSEMBL"
    symbol <- "SYMBOL"
    txpre <- "ENSRNOT"
    gxpre <- "ENSRNOG"
  } else if (organism == "Celegans") {
    package <- "org.Ce.eg.db"
    keytype <- "WORMBASE"
    symbol <- "SYMBOL"
    txpre <- ""
    gxpre <- "WBGene"
  } else if (organism == "Dmelanogaster") {
    package <- "org.Dm.eg.db"
    keytype <- "FLYBASE"
    symbol <- "SYMBOL"
    txpre <- "FBtr"
    gxpre <- "FBgn"
  } else if (organism == "Scerevisiae") {
    package <- "org.Sc.sgd.db"
    keytype <- "GENENAME"
    symbol <- "GENENAME"
    txpre <- ""
    gxpre <- ""
  } else {
    stop("Unable to find annotation support for",organism,"...patches welcome!")
  }

  org <- list(package=package, 
              txpre=txpre, 
              gxpre=gxpre, 
              keytype=keytype,
              symbol=symbol)
  return(org)

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

  if (type == "ErccDbLite") {
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
#' 
#' @return the name of the annotation package, or NULL if uncertain how to do it
#' 
#' @export
#'
createAnnotationPackage <- function(fastaFile, ...) { # {{{

  type <- getAnnotationType(fastaFile) 
  if (type == "EnsDbLite") { 
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
