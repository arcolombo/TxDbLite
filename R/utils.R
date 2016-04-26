#' utility functions
#'
utils <- NULL


#' @describeIn utils
#' 
#' needed from Ens83 onwards :-/ and maybe will need to generalize as well
#' 
#' @param x   a string to try and pop the ENSEMBL identifier off of
#'
#' @return    if everything works out properly, the ENS* identifier (without .N)
#'
#' @export
idStub <- function(x) {
  if (grepl("^ENS.*\\.", x)) x <- strsplit(x, "\\.")[[1]][1]
  if (grepl(" ", x)) x <- strsplit(x, " ")[[1]][1]
  return(x)
}


#' @describeIn utils
#' 
#' get the "stub" of a FASTA filename (no .fa, no .fasta, no .gz)
#' 
#' @export
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
colAsVector <- function(x, y) t(x[, y, drop=FALSE])[1,]

#' @describeIn utils 
#' 
#' supported organism abbreviations
#'
#' @export
getSupportedAbbreviations <- function(how=c("abbreviation","reactome")) { # {{{
  how <- match.arg(how) 
  if (!exists("supportedOrganismsForTxDbLite")) 
    data(supportedOrganismsForTxDbLite, package="TxDbLite")
  colAsVector(supportedOrganismsForTxDbLite, how)
} # }}}

#' @describeIn utils
#' 
#' Get the abbreviation for an organism (among those with which we're familiar)
#' 
#' @export
getOrganismAbbreviation <- function(organism, how=c("abbreviation","reactome")){# {{{

  how <- match.arg(how)
  org <- sub(" ", "_", sub("\\.", "_", organism))
  abbr <- getSupportedAbbreviations(how=how)
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
    organism <- sub("\\.", "_", # try to be robust
                    sub("Mmusculus", "Mus_musculus", 
                        sub("Hsapiens", "Homo_sapiens",
                          sub("Dmelanogaster","Drosophila_melanogaster", 
                              organism))))

    organism <- getOrganismAbbreviation(organism)
    genomeVersion <- tokens[2]
    if (length(tokens) > 3) {
      version <- tokens[3]
      what <- tokens[4]
      if(grepl("_", what)) {
         what <- gsub("_",".",what)
      }
    
    } else { 
      version <- tokens[2]
      what <- tokens[3] 
      if(grepl("_",what)) { 
        what <- gsub("_", ".", what)
      }
    }
    return(gsub("_", "", paste(type, organism, version, sep=".")))
  } else {
    return(NULL)
  }

} # }}}
